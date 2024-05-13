// This is an advanced implementation of the algorithm described in the following paper:
//   J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time.
//     Robotics: Science and Systems Conference (RSS). Berkeley, CA, July 2014.

// Modifier: Tong Qin               qintonguav@gmail.com
// 	         Shaozu Cao 		    saozu.cao@connect.ust.hk

// Copyright 2013, Ji Zhang, Carnegie Mellon University
// Further contributions copyright (c) 2016, Southwest Research Institute
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <cmath>
#include <vector>
#include <string>
#include "aloam_velodyne/common.h"
#include "aloam_velodyne/tic_toc.h"
#include <nav_msgs/Odometry.h>
#include <opencv/cv.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>

using std::atan2;
using std::cos;
using std::sin;

// 扫描周期, velodyne频率10Hz，周期0.1s
const double scanPeriod = 0.1;

// 弃用前systemDelay帧初始数据
const int systemDelay = 0;

// systemInitCount用于计数过了多少帧
// 超过systemDelay后，systemInited为true即初始化完成
int systemInitCount = 0;
bool systemInited = false;
int N_SCANS = 0;

// 点云曲率, 400000为一帧点云中点的最大数量
float cloudCurvature[400000];

// 曲率点对应的序号
int cloudSortInd[400000];

// 点是否筛选过标志：0-未筛选过，1-筛选过
int cloudNeighborPicked[400000];

// 点分类标号:2-代表曲率很大，1-代表曲率比较大, 0-曲率比较小, -1-代表曲率很小(其中1包含了2, 0包含了1, 0和1构成了点云全部的点)
int cloudLabel[400000];

// 两点曲率比较
bool comp(int i, int j) { return (cloudCurvature[i] < cloudCurvature[j]); }

ros::Publisher pubLaserCloud;
ros::Publisher pubCornerPointsSharp;
ros::Publisher pubCornerPointsLessSharp;
ros::Publisher pubSurfPointsFlat;
ros::Publisher pubSurfPointsLessFlat;
ros::Publisher pubRemovePoints;

// ros形式的一线扫描
std::vector<ros::Publisher> pubEachScan;

bool PUB_EACH_LINE = false;

// 根据距离去除过近的点，距离的参数
double MINIMUM_RANGE = 0.1;

// 去除半径thres内的点
template <typename PointT>
void removeClosedPointCloud(const pcl::PointCloud<PointT> &cloud_in,
                            pcl::PointCloud<PointT> &cloud_out, float thres)
{
    if (&cloud_in != &cloud_out)
    {
        cloud_out.header = cloud_in.header;
        cloud_out.points.resize(cloud_in.points.size());
    }

    size_t j = 0;

    for (size_t i = 0; i < cloud_in.points.size(); ++i)
    {
        if (cloud_in.points[i].x * cloud_in.points[i].x + cloud_in.points[i].y * cloud_in.points[i].y + cloud_in.points[i].z * cloud_in.points[i].z < thres * thres)
            continue;
        cloud_out.points[j] = cloud_in.points[i];
        j++;
    }
    if (j != cloud_in.points.size())
    {
        cloud_out.points.resize(j);
    }

    cloud_out.height = 1;
    cloud_out.width = static_cast<uint32_t>(j);
    cloud_out.is_dense = true;
}

// 订阅点云句柄
void laserCloudHandler(const sensor_msgs::PointCloud2ConstPtr &laserCloudMsg)
{
    // 剔除前systemDelay帧
    if (!systemInited)
    {
        systemInitCount++;
        if (systemInitCount >= systemDelay)
        {
            systemInited = true;
        }
        else
            return;
    }

    TicToc t_whole;
    TicToc t_prepare;
    // 记录每个scan有曲率的点的开始和结束索引
    std::vector<int> scanStartInd(N_SCANS, 0);
    std::vector<int> scanEndInd(N_SCANS, 0);

    pcl::PointCloud<pcl::PointXYZ> laserCloudIn;
    // 把ros包的点云转化为pcl形式
    pcl::fromROSMsg(*laserCloudMsg, laserCloudIn);
    std::vector<int> indices;

    // 这个函数目的是去除无效点，第一个参数是输入，第二个参数是输出，第三个参数是输出的点在输入里的位置
    // cloud_out.points[i] = cloud_in.points[indices[i]]
    pcl::removeNaNFromPointCloud(laserCloudIn, laserCloudIn, indices);
    removeClosedPointCloud(laserCloudIn, laserCloudIn, MINIMUM_RANGE);

    int cloudSize = laserCloudIn.points.size();
    // 每次扫描是一条线，看作者的数据集激光x向前，y向左，那么下面就是线一端到另一端
    // 计算旋转角时取负号是因为velodyne是顺时针旋转
    float startOri = -atan2(laserCloudIn.points[0].y, laserCloudIn.points[0].x);
    float endOri = -atan2(laserCloudIn.points[cloudSize - 1].y,
                          laserCloudIn.points[cloudSize - 1].x) +
                   2 * M_PI;

    // 激光间距收束到 pi < endOri - startOri < 3 * pi
    if (endOri - startOri > 3 * M_PI)
    {
        endOri -= 2 * M_PI;
    }
    else if (endOri - startOri < M_PI)
    {
        endOri += 2 * M_PI;
    }
    // printf("end Ori %f\n", endOri);

    // 过半记录标志
    bool halfPassed = false;
    // 记录总点数
    int count = cloudSize;
    PointType point;
    std::vector<pcl::PointCloud<PointType>> laserCloudScans(N_SCANS);
    for (int i = 0; i < cloudSize; i++)
    {
        point.x = laserCloudIn.points[i].x;
        point.y = laserCloudIn.points[i].y;
        point.z = laserCloudIn.points[i].z;

        // 求仰角atan输出为-pi/2到pi/2
        float angle = atan(point.z / sqrt(point.x * point.x + point.y * point.y)) * 180 / M_PI;
        // printf("angle = %.2f\n", angle);
        int scanID = 0;

        if (N_SCANS == 16)
        {
            scanID = int((angle + 15) / 2 + 0.5);     // 属于第几根线，+0.5再int，相当于四舍五入操作，垂直分辨率为2°，范围为-15°~15°
            if (scanID > (N_SCANS - 1) || scanID < 0) // 计算得到的scanID不合法
            {
                count--;
                continue;
            }
        }
        else if (N_SCANS == 32)
        {
            scanID = int((angle + 92.0 / 3.0) * 3.0 / 4.0);
            if (scanID > (N_SCANS - 1) || scanID < 0)
            {
                count--;
                continue;
            }
        }
        else if (N_SCANS == 64)
        {
            if (angle >= -8.83)
                scanID = int((2 - angle) * 3.0 + 0.5);
            else
                scanID = N_SCANS / 2 + int((-8.83 - angle) * 2.0 + 0.5);

            // use [0 50]  > 50 remove outlies
            if (angle > 2 || angle < -24.33 || scanID > 50 || scanID < 0)
            {
                count--;
                continue;
            }
        }
        else if (N_SCANS == 128)
        {
            scanID = int((angle + 22.15) / 0.35 + 0.5);
            if (scanID > (N_SCANS - 1) || scanID < 0)
            {
                count--;
                continue;
            }
        }
        else
        {
            printf("wrong scan number\n");
            ROS_BREAK();
        }
        // printf("angle %f scanID %d \n", anSubscgle, scanID);

        // 当前旋转角
        float ori = -atan2(point.y, point.x);
        // 根据扫描线是否旋转过半选择与起始位置还是终止位置进行差值计算，从而进行补偿
        if (!halfPassed)
        {
            // 确保-pi/2 < ori - startOri < 3*pi/2
            if (ori < startOri - M_PI / 2)
            {
                ori += 2 * M_PI;
            }
            else if (ori > startOri + M_PI * 3 / 2)
            {
                ori -= 2 * M_PI;
            }
            // 旋转过半
            if (ori - startOri > M_PI)
            {
                halfPassed = true;
            }
        }
        else
        {
            // 确保-pi/2 < endOri - ori < 3*pi/2
            ori += 2 * M_PI;
            if (ori < endOri - M_PI * 3 / 2)
            {
                ori += 2 * M_PI;
            }
            else if (ori > endOri + M_PI / 2)
            {
                ori -= 2 * M_PI;
            }
        }
        // 计算每一个点由于旋转导致的时间滞后relTime
        float relTime = (ori - startOri) / (endOri - startOri);
        // 用point.intensity记录当前线束ID和本线扫描时间之和
        point.intensity = scanID + scanPeriod * relTime;
        // 按线分类保存
        laserCloudScans[scanID].push_back(point);
    }

    cloudSize = count;
    printf("points size %d \n", cloudSize);

    // 把所有线保存在laserCloud一个数据集合里，把每条线的第五个和倒数第五个位置反馈给scanStartInd和scanEndInd
    pcl::PointCloud<PointType>::Ptr laserCloud(new pcl::PointCloud<PointType>());
    // 去除前后5个点
    for (int i = 0; i < N_SCANS; i++)
    {
        scanStartInd[i] = laserCloud->size() + 5;
        *laserCloud += laserCloudScans[i];
        scanEndInd[i] = laserCloud->size() - 6;
    }

    printf("prepare time %f \n", t_prepare.toc());
    // 求曲率，这里选取一条线上的十个点
    for (int i = 5; i < cloudSize - 5; i++)
    {
        float diffX = laserCloud->points[i - 5].x + laserCloud->points[i - 4].x + laserCloud->points[i - 3].x + laserCloud->points[i - 2].x + laserCloud->points[i - 1].x - 10 * laserCloud->points[i].x + laserCloud->points[i + 1].x + laserCloud->points[i + 2].x + laserCloud->points[i + 3].x + laserCloud->points[i + 4].x + laserCloud->points[i + 5].x;
        float diffY = laserCloud->points[i - 5].y + laserCloud->points[i - 4].y + laserCloud->points[i - 3].y + laserCloud->points[i - 2].y + laserCloud->points[i - 1].y - 10 * laserCloud->points[i].y + laserCloud->points[i + 1].y + laserCloud->points[i + 2].y + laserCloud->points[i + 3].y + laserCloud->points[i + 4].y + laserCloud->points[i + 5].y;
        float diffZ = laserCloud->points[i - 5].z + laserCloud->points[i - 4].z + laserCloud->points[i - 3].z + laserCloud->points[i - 2].z + laserCloud->points[i - 1].z - 10 * laserCloud->points[i].z + laserCloud->points[i + 1].z + laserCloud->points[i + 2].z + laserCloud->points[i + 3].z + laserCloud->points[i + 4].z + laserCloud->points[i + 5].z;

        cloudCurvature[i] = diffX * diffX + diffY * diffY + diffZ * diffZ;
        cloudSortInd[i] = i;        // 排序时用到的索引
        cloudNeighborPicked[i] = 0; // 是否筛选过
        cloudLabel[i] = 0;          // 标记点的类型
    }

    TicToc t_pts;

    // 角点，降采样角点，面点，降采样面点
    pcl::PointCloud<PointType> cornerPointsSharp;
    pcl::PointCloud<PointType> cornerPointsLessSharp;
    pcl::PointCloud<PointType> surfPointsFlat;
    pcl::PointCloud<PointType> surfPointsLessFlat;

    float t_q_sort = 0;
    for (int i = 0; i < N_SCANS; i++)
    {
        // 以6根扫描线束为界限来计算,点数小于6就认为是边界部分，会跳过本次循环
        if (scanEndInd[i] - scanStartInd[i] < 6)
            continue;
        pcl::PointCloud<PointType>::Ptr surfPointsLessFlatScan(new pcl::PointCloud<PointType>);
        for (int j = 0; j < 6; j++)
        {
            // 六等份起点：sp = scanStartInd + (scanEndInd - scanStartInd)*j/6
            // 六等份终点：ep = scanStartInd - 1 + (scanEndInd - scanStartInd)*(j+1)/6
            int sp = scanStartInd[i] + (scanEndInd[i] - scanStartInd[i]) * j / 6;
            int ep = scanStartInd[i] + (scanEndInd[i] - scanStartInd[i]) * (j + 1) / 6 - 1;

            TicToc t_tmp;
            std::sort(cloudSortInd + sp, cloudSortInd + ep + 1, comp);
            t_q_sort += t_tmp.toc();
            // 挑选每个分段的曲率很大和比较大的点
            int largestPickedNum = 0;
            for (int k = ep; k >= sp; k--)
            {
                int ind = cloudSortInd[k];
                // 如果筛选标志为0，并且曲率较大
                if (cloudNeighborPicked[ind] == 0 &&
                    cloudCurvature[ind] > 0.1)
                {
                    // 则曲率大的点记数一下
                    largestPickedNum++;
                    if (largestPickedNum <= 2)
                    {
                        cloudLabel[ind] = 2;
                        cornerPointsSharp.push_back(laserCloud->points[ind]);
                        cornerPointsLessSharp.push_back(laserCloud->points[ind]);
                    }
                    // 保存20个角点
                    else if (largestPickedNum <= 20)
                    {
                        cloudLabel[ind] = 1;
                        cornerPointsLessSharp.push_back(laserCloud->points[ind]);
                    }
                    else
                    {
                        break;
                    }

                    cloudNeighborPicked[ind] = 1;
                    // 将曲率比较大的点的前后各5个连续距离比较近的点筛选出去，防止特征点聚集，使得特征点在每个方向上尽量分布均匀
                    for (int l = 1; l <= 5; l++)
                    {
                        float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l - 1].x;
                        float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l - 1].y;
                        float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l - 1].z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                        {
                            break;
                        }

                        cloudNeighborPicked[ind + l] = 1;
                    }
                    for (int l = -1; l >= -5; l--)
                    {
                        float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l + 1].x;
                        float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l + 1].y;
                        float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l + 1].z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                        {
                            break;
                        }

                        cloudNeighborPicked[ind + l] = 1;
                    }
                }
            }

            // 挑选每个分段的曲率很小比较小的点
            int smallestPickedNum = 0;
            // 升序找到最小的点，来作为平面点
            for (int k = sp; k <= ep; k++)
            {
                int ind = cloudSortInd[k];

                if (cloudNeighborPicked[ind] == 0 &&
                    cloudCurvature[ind] < 0.1)
                {

                    cloudLabel[ind] = -1;
                    surfPointsFlat.push_back(laserCloud->points[ind]);
                    // 只选最小的四个，剩下的Label==0,就都是曲率比较小的
                    smallestPickedNum++;
                    if (smallestPickedNum >= 4)
                    {
                        break;
                    }

                    cloudNeighborPicked[ind] = 1;
                    for (int l = 1; l <= 5; l++)
                    {
                        float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l - 1].x;
                        float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l - 1].y;
                        float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l - 1].z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                        {
                            break;
                        }

                        cloudNeighborPicked[ind + l] = 1;
                    }
                    for (int l = -1; l >= -5; l--)
                    {
                        float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l + 1].x;
                        float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l + 1].y;
                        float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l + 1].z;
                        if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
                        {
                            break;
                        }

                        cloudNeighborPicked[ind + l] = 1;
                    }
                }
            }

            // 将剩余的点（包括之前被排除的点）全部归入平面点中less flat类别中
            for (int k = sp; k <= ep; k++)
            {
                if (cloudLabel[k] <= 0)
                {
                    surfPointsLessFlatScan->push_back(laserCloud->points[k]);
                }
            }
        }
        // 由于less flat点最多，对每个分段less flat的点进行体素栅格滤波（下采样）
        pcl::PointCloud<PointType> surfPointsLessFlatScanDS;
        pcl::VoxelGrid<PointType> downSizeFilter;
        downSizeFilter.setInputCloud(surfPointsLessFlatScan);
        downSizeFilter.setLeafSize(0.2, 0.2, 0.2);
        downSizeFilter.filter(surfPointsLessFlatScanDS);

        surfPointsLessFlat += surfPointsLessFlatScanDS;
    }
    printf("sort q time %f \n", t_q_sort);
    printf("seperate points time %f \n", t_pts.toc());
    // 发布信息准备工作，主要是整理后的点云、角点面点及其降采样结果
    sensor_msgs::PointCloud2 laserCloudOutMsg;
    pcl::toROSMsg(*laserCloud, laserCloudOutMsg);
    laserCloudOutMsg.header.stamp = laserCloudMsg->header.stamp;
    laserCloudOutMsg.header.frame_id = "/camera_init";
    pubLaserCloud.publish(laserCloudOutMsg);

    sensor_msgs::PointCloud2 cornerPointsSharpMsg;
    pcl::toROSMsg(cornerPointsSharp, cornerPointsSharpMsg);
    cornerPointsSharpMsg.header.stamp = laserCloudMsg->header.stamp;
    cornerPointsSharpMsg.header.frame_id = "/camera_init";
    pubCornerPointsSharp.publish(cornerPointsSharpMsg);

    sensor_msgs::PointCloud2 cornerPointsLessSharpMsg;
    pcl::toROSMsg(cornerPointsLessSharp, cornerPointsLessSharpMsg);
    cornerPointsLessSharpMsg.header.stamp = laserCloudMsg->header.stamp;
    cornerPointsLessSharpMsg.header.frame_id = "/camera_init";
    pubCornerPointsLessSharp.publish(cornerPointsLessSharpMsg);

    sensor_msgs::PointCloud2 surfPointsFlat2;
    pcl::toROSMsg(surfPointsFlat, surfPointsFlat2);
    surfPointsFlat2.header.stamp = laserCloudMsg->header.stamp;
    surfPointsFlat2.header.frame_id = "/camera_init";
    pubSurfPointsFlat.publish(surfPointsFlat2);

    sensor_msgs::PointCloud2 surfPointsLessFlat2;
    pcl::toROSMsg(surfPointsLessFlat, surfPointsLessFlat2);
    surfPointsLessFlat2.header.stamp = laserCloudMsg->header.stamp;
    surfPointsLessFlat2.header.frame_id = "/camera_init";
    pubSurfPointsLessFlat.publish(surfPointsLessFlat2);

    // pub each scam
    if (PUB_EACH_LINE)
    {
        for (int i = 0; i < N_SCANS; i++)
        {
            sensor_msgs::PointCloud2 scanMsg;
            pcl::toROSMsg(laserCloudScans[i], scanMsg);
            scanMsg.header.stamp = laserCloudMsg->header.stamp;
            scanMsg.header.frame_id = "/camera_init";
            pubEachScan[i].publish(scanMsg);
        }
    }

    printf("scan registration time %f ms *************\n", t_whole.toc());
    if (t_whole.toc() > 100)
        ROS_WARN("scan registration process over 100ms");
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "scanRegistration");
    ros::NodeHandle nh;
    // 参数，线数
    nh.param<int>("scan_line", N_SCANS, 16);
    nh.param<double>("minimum_range", MINIMUM_RANGE, 0.1);

    printf("scan line number %d \n", N_SCANS);

    if (N_SCANS != 16 && N_SCANS != 32 && N_SCANS != 64 && N_SCANS != 128)
    {
        printf("only support velodyne with 16, 32, 64 or 128 scan line!");
        return 0;
    }
    // 接受硬件传来的点云数据
    ros::Subscriber subLaserCloud = nh.subscribe<sensor_msgs::PointCloud2>("/ouster/points", 100, laserCloudHandler);
    
    // 发布去除无效点、噪声点，并进行扫描线号分开存储的点云数据
    pubLaserCloud = nh.advertise<sensor_msgs::PointCloud2>("/velodyne_cloud_2", 100);
    // 发布提取出的特征点，包括以下:角点，降采样角点，面点，降采样面点
    pubCornerPointsSharp = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_sharp", 100);
    pubCornerPointsLessSharp = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_less_sharp", 100);
    pubSurfPointsFlat = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_flat", 100);
    pubSurfPointsLessFlat = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_less_flat", 100);
    // 发布被去除掉的点
    pubRemovePoints = nh.advertise<sensor_msgs::PointCloud2>("/laser_remove_points", 100);
    // 发布每行扫描线
    if (PUB_EACH_LINE)
    {
        for (int i = 0; i < N_SCANS; i++)
        {
            ros::Publisher tmp = nh.advertise<sensor_msgs::PointCloud2>("/laser_scanid_" + std::to_string(i), 100);
            pubEachScan.push_back(tmp);
        }
    }
    ros::spin();

    return 0;
}
