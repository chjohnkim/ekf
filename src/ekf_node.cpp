#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/Range.h>
#include <nav_msgs/Odometry.h>
#include "visualization_msgs/Marker.h"
#include <Eigen/Eigen>
#define gravity 9.81

using namespace std;
using namespace Eigen;
ros::Publisher ekf_odom_pub;
ros::Publisher ref_odom_pub;
ros::Publisher opflow_vel_pub;
ros::Publisher vicon_vel_pub;
ros::Publisher ekf_vel_pub;
MatrixXd Q = MatrixXd::Identity(12, 12);
MatrixXd Rt = MatrixXd::Identity(9,9);
MatrixXd Rt2 = MatrixXd::Identity(4,4);

VectorXd mean = Eigen::VectorXd::Zero(15);
MatrixXd var  = Eigen::MatrixXd::Zero(15,15);
VectorXd state_measured = Eigen::VectorXd::Zero(15);
VectorXd latest_mean = Eigen::VectorXd::Zero(15);
MatrixXd latest_var  = Eigen::MatrixXd::Zero(15,15);

MatrixXd A = MatrixXd::Zero(15,15);
MatrixXd U = MatrixXd::Zero(15,12);
MatrixXd F = MatrixXd::Zero(15,15);
MatrixXd V = MatrixXd::Zero(15,12);
MatrixXd I3 = MatrixXd::Identity(3,3);
MatrixXd I15 = MatrixXd::Identity(15,15);

VectorXd g = Eigen::VectorXd::Zero(3);;
Eigen::Vector3d last_rpy;


double tag_odom_time, opflow_odom_time, imu_time, latest_time;
bool initialized = false; 
double tag_count = 0, last_tag_count = -1;


void visualizeVelocity(Eigen::Vector3d position, Eigen::Vector3d velocity,
                       int id, Eigen::Vector3d color, ros::Publisher pub_vel) {
    double scale = 10;
    visualization_msgs::Marker m;
    m.header.frame_id = "world";
    m.id = id;
    m.type = visualization_msgs::Marker::ARROW;
    m.action = visualization_msgs::Marker::MODIFY;
    m.scale.x = 0.1;
    m.scale.y = 0.3;
    m.scale.z = 0;
    m.pose.position.x = 0;
    m.pose.position.y = 0;
    m.pose.position.z = 0;
    m.pose.orientation.w = 1;
    m.pose.orientation.x = 0;
    m.pose.orientation.y = 0;
    m.pose.orientation.z = 0;
    m.color.a = 1.0;
    m.color.r = color.x();
    m.color.g = color.y();
    m.color.b = color.z();
    m.points.clear();
    geometry_msgs::Point point;
    point.x = position.x();
    point.y = position.y();
    point.z = position.z();
    m.points.push_back(point);
    point.x = position.x() + velocity.x() * scale;
    point.y = position.y() + velocity.y() * scale;
    point.z = position.z() + velocity.z() * scale;
    m.points.push_back(point);
    pub_vel.publish(m);
}


Eigen::MatrixXd Ginv_dotproduct_w(double phi, double theta, double psi, Vector3d w)
{
    Eigen::MatrixXd Ginv_dot_w = Eigen::MatrixXd::Zero(3,3);
    Ginv_dot_w << 0,                                                                w(2)*cos(theta) - w(0)*sin(theta),                      0,
                -(w(2)*cos(theta) - w(0)*sin(theta))/pow(cos(phi),2),              (sin(phi)*(w(0)*cos(theta) + w(2)*sin(theta)))/cos(phi), 0,
                 (sin(phi)*(w(2)*cos(theta) - w(0)*sin(theta)))/pow(cos(phi),2),  -(w(0)*cos(theta) + w(2)*sin(theta))/cos(phi),            0;

    return Ginv_dot_w;
}
     
Eigen::MatrixXd R_dotproduct_acc(double phi, double theta, double psi, Vector3d acc)
{
    Eigen::MatrixXd R_dot_acc = Eigen::MatrixXd::Zero(3,3);
    R_dot_acc  << sin(psi)*(acc(1)*sin(phi) + acc(2)*cos(phi)*cos(theta) - acc(0)*cos(phi)*sin(theta)), acc(2)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) - acc(0)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)),  -acc(0)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - acc(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) - acc(1)*cos(phi)*cos(psi),
                 -cos(psi)*(acc(1)*sin(phi) + acc(2)*cos(phi)*cos(theta) - acc(0)*cos(phi)*sin(theta)), acc(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - acc(0)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)),   acc(0)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + acc(2)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) - acc(1)*cos(phi)*sin(psi),
                  acc(1)*cos(phi) - acc(2)*cos(theta)*sin(phi) + acc(0)*sin(phi)*sin(theta),           -cos(phi)*(acc(0)*cos(theta) + acc(2)*sin(theta)),                                                                            0;

    return R_dot_acc;
}
 
Eigen::Matrix3d euler_to_R(double phi, double theta, double psi)
{
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3,3);
    R(0,0) =  cos(psi)*cos(theta) - sin(phi)*sin(theta)*sin(psi);
    R(0,1) = -cos(phi)*sin(psi);
    R(0,2) =  cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi);
    R(1,0) =  cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta);
    R(1,1) =  cos(phi)*cos(psi);
    R(1,2) =  sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi);
    R(2,0) = -cos(phi)*sin(theta);  
    R(2,1) =  sin(phi);
    R(2,2) =  cos(phi)*cos(theta); 
    return R;
}

MatrixXd euler_to_G(double phi, double theta, double psi)
{
    MatrixXd G = MatrixXd::Zero(3, 3);
    G <<    cos(theta),    0,     -cos(phi)*sin(theta),
            0,             1,           sin(phi),
            sin(theta),    0,      cos(phi)*cos(theta);
    return G;
}

MatrixXd euler_to_G_inverse(double phi, double theta, double psi)
{
    MatrixXd Ginv = MatrixXd::Zero(3, 3);
    Ginv<<  cos(theta),                       0,               sin(theta),
            (sin(phi)*sin(theta))/cos(phi),    1,              -(cos(theta)*sin(phi))/cos(phi),
            -sin(theta)/cos(phi),              0,               cos(theta)/cos(phi);
    return Ginv;
}

Eigen::Vector3d R_to_euler(const Eigen::Matrix3d& R)
{
  Eigen::Vector3d euler;
  double r = asin(R(2,1));
  double p = atan2(-R(2,0)/cos(r), R(2,2)/cos(r));
  double y = atan2(-R(0,1)/cos(r), R(1,1)/cos(r));
  euler(0) = r;
  euler(1) = p;
  euler(2) = y;

  return euler;
}


void viconCallback(const nav_msgs::Odometry::ConstPtr &vicon_msg) {
    Eigen::Vector3d vicon_vel, vicon_pos;
    vicon_pos.x() = vicon_msg->pose.pose.position.x;
    vicon_pos.y() = vicon_msg->pose.pose.position.y;
    vicon_pos.z() = vicon_msg->pose.pose.position.z;
    vicon_vel.x() = vicon_msg->twist.twist.linear.x;
    vicon_vel.y() = vicon_msg->twist.twist.linear.y;
    vicon_vel.z() = vicon_msg->twist.twist.linear.z;
    visualizeVelocity(vicon_pos, vicon_vel, 1, Eigen::Vector3d(0,0,1), vicon_vel_pub);
}

void imu_callback(const sensor_msgs::Imu::ConstPtr &msg)
{
    //your code for propagation
    imu_time = msg->header.stamp.toSec();
    Eigen::Vector3d ang_vel_imu, lin_acc_imu;
    ang_vel_imu.x() = msg->angular_velocity.x;
    ang_vel_imu.y() = msg->angular_velocity.y;
    ang_vel_imu.z() = msg->angular_velocity.z;
    lin_acc_imu.x() = msg->linear_acceleration.x;
    lin_acc_imu.y() = msg->linear_acceleration.y;
    lin_acc_imu.z() = msg->linear_acceleration.z;

    if (initialized == false)
    {
        cout << "Initializing......" << endl;
        return;
    }
    
    //PROPAGATION
    double dt = imu_time-latest_time;
    Vector3d X1, X2, X3, X4, X5;
    X1 = latest_mean.segment<3>(0); // POSITION
    X2 = latest_mean.segment<3>(3); // ORIENTATION
    X3 = latest_mean.segment<3>(6); // VELOCITY
    X4 = latest_mean.segment<3>(9); // NOISE (?)
    X5 = latest_mean.segment<3>(12); // NOISE (?)
    
    double phi = X2(0), theta = X2(1), psi = X2(2);
    Vector3d w   = ang_vel_imu - X4;
    Vector3d acc = lin_acc_imu - X5;
    
    A.block<3, 3>(0, 6) =  I3;
    A.block<3, 3>(3, 3) = Ginv_dotproduct_w(phi,theta,psi,w); 
    A.block<3, 3>(3, 9) = -euler_to_G(phi,theta,psi).inverse(); 
   	A.block<3, 3>(6, 3) = R_dotproduct_acc(phi,theta,psi,acc); 
   	A.block<3, 3>(6,12) = -euler_to_R(phi,theta,psi);

    U.block<3, 3>(3, 0) = -euler_to_G(phi,theta,psi).inverse();
    U.block<3, 3>(6, 3) = -euler_to_R(phi,theta,psi);
    U.block<3, 3>(9, 6) = I3;
    U.block<3, 3>(12,9) = I3;
    F = I15 + dt * A;
    V = dt * U;

    VectorXd xdot = Eigen::VectorXd::Zero(15);
    xdot.segment(0, 3) = X3;
    xdot.segment(3, 3) = euler_to_G_inverse(phi,theta,psi) * w;
    xdot.segment(6, 3) = euler_to_R(phi,theta,psi) * acc + g;

    latest_mean = latest_mean + dt * xdot;
    latest_var = F*latest_var*F.transpose() + V*Q*V.transpose(); 
    latest_time = imu_time;
    //cout << latest_mean(3) << "   " << latest_mean(3) << "   " << latest_mean(3) << endl;
    //Publish propgation odometry
    Quaterniond Q_ekf;
    Q_ekf = euler_to_R(latest_mean(3),latest_mean(4),latest_mean(5));
    nav_msgs::Odometry ekf_odom;
    ekf_odom.header.stamp = msg->header.stamp;
    ekf_odom.header.frame_id = "world";
    ekf_odom.pose.pose.position.x = latest_mean(0);
    ekf_odom.pose.pose.position.y = latest_mean(1);
    ekf_odom.pose.pose.position.z = latest_mean(2);
    ekf_odom.pose.pose.orientation.w = Q_ekf.w();
    ekf_odom.pose.pose.orientation.x = Q_ekf.x();
    ekf_odom.pose.pose.orientation.y = Q_ekf.y();
    ekf_odom.pose.pose.orientation.z = Q_ekf.z();
    ekf_odom.header.stamp = msg->header.stamp;
    ekf_odom.header.frame_id = "world";
    ekf_odom.twist.twist.linear.x = latest_mean(6);
    ekf_odom.twist.twist.linear.y = latest_mean(7);
    ekf_odom.twist.twist.linear.z = latest_mean(8);
    ekf_odom_pub.publish(ekf_odom);

    visualizeVelocity(Eigen::Vector3d(latest_mean(0), latest_mean(1), latest_mean(2)), Eigen::Vector3d(latest_mean(6), latest_mean(7), latest_mean(8)), 2, Eigen::Vector3d(1,0,0), ekf_vel_pub);

}

//Rotation from the camera frame to the IMU frame
Eigen::Matrix3d Rcam;
void odom_callback_tag(const nav_msgs::Odometry::ConstPtr &msg)
{
    tag_odom_time = msg->header.stamp.toSec();
    Eigen::Quaterniond q_cw;
    Vector3d T_cw, T_ic, T_wi, euler_wi, T_offset;
    Matrix3d R_cw, R_ic, R_wi;
    T_cw.x() = msg->pose.pose.position.x;
    T_cw.y() = msg->pose.pose.position.y;
    T_cw.z() = msg->pose.pose.position.z;
    q_cw = Eigen::Quaterniond(msg->pose.pose.orientation.w,
                           msg->pose.pose.orientation.x,
                           msg->pose.pose.orientation.y,
                           msg->pose.pose.orientation.z);
    
    MatrixXd R_vm = MatrixXd::Zero(3, 3);
    R_vm << 0,    1,     0,
            1,    0,     0,
            0,    0,    -1;
    
    // Camera IMU pose for Project 2 Phase 2 Part 2 & 3
    T_ic << -0.1, 0.0, -0.03;
    R_ic = Quaterniond(0, -0.7071068, 0.7071068, 0).toRotationMatrix();
    
    /*
    R_cw = q_cw.toRotationMatrix();
    R_wi = Quaterniond(0, 0.7071068, 0.7071068, 0).toRotationMatrix()*((R_ic*R_cw).transpose());
    T_wi = Quaterniond(0, 0.7071068, 0.7071068, 0).toRotationMatrix()*(-R_cw.transpose() *(R_ic.transpose() * T_ic + T_cw)); 
    euler_wi = R_to_euler(R_wi);
    */
    
    T_offset << -0.5, -0.5, 0;

    R_cw = q_cw.toRotationMatrix();
    R_wi = (R_ic*R_cw).transpose();
    R_wi = R_vm*R_wi;
    T_wi = -R_cw.transpose() *(R_ic.transpose() * T_ic + T_cw);
    T_wi = R_vm*T_wi + T_offset; 
    euler_wi = R_to_euler(R_wi);
    

    
    if (initialized == false)
    {
        g << 0, 0, -gravity;
        latest_mean.segment<3>(0) << T_wi;
        latest_mean.segment<3>(3) << euler_wi;
        latest_time = tag_odom_time;
        initialized = true;
        cout << "Initialized!" << endl;
        return;
    }
      
    state_measured.segment<3>(0) << T_wi;
    state_measured.segment<3>(3) << euler_wi;
    tag_count = tag_count + 1;

    //Publish tag odometry 
    Quaterniond Q_ref;
    Q_ref = R_wi;
    nav_msgs::Odometry ref_odom;
    ref_odom.header.stamp = msg->header.stamp;
    ref_odom.header.frame_id = "world";
    ref_odom.pose.pose.position.x = T_wi(0);
    ref_odom.pose.pose.position.y = T_wi(1);
    ref_odom.pose.pose.position.z = T_wi(2);
    ref_odom.pose.pose.orientation.w = Q_ref.w();
    ref_odom.pose.pose.orientation.x = Q_ref.x();
    ref_odom.pose.pose.orientation.y = Q_ref.y();
    ref_odom.pose.pose.orientation.z = Q_ref.z();
    ref_odom_pub.publish(ref_odom);
    
    /*
    //UPDATE FOR TAG ODOM ONLY
    MatrixXd C = MatrixXd::Identity(6,15);
    MatrixXd W = MatrixXd::Identity(6,6);
    MatrixXd K = latest_var*C.transpose()*(C*latest_var*C.transpose() + W*(Rt.block<6, 6>(0, 0))*W.transpose()).inverse();

    
    for(int i=3;i<6;i++){
        if(abs(last_rpy(i-3)-state_measured(i))>M_PI){
            if(state_measured(i)<0) state_measured(i) = 2 * M_PI + state_measured(i);
            else state_measured(i) = -2* M_PI + state_measured(i);   
        }         
    }
    last_rpy(0) = state_measured(3);
    last_rpy(1) = state_measured(4);
    last_rpy(2) = state_measured(5);
    
    latest_mean = latest_mean + K*(C*state_measured - C*latest_mean);
    latest_var = latest_var - K*C*latest_var;
    */
}


void odom_callback_opflow(const nav_msgs::Odometry::ConstPtr &msg)
{
    opflow_odom_time = msg->header.stamp.toSec();
    Matrix3d R_wc;
    Vector3d V_bc, V_wc;

    
    V_bc.x() = msg->twist.twist.linear.x;
    V_bc.y() = msg->twist.twist.linear.y;
    V_bc.z() = msg->twist.twist.linear.z;
    double height = msg->pose.pose.position.z; 
    R_wc = Quaterniond(0, 0.7071068, 0.7071068, 0).toRotationMatrix();
    //V_wc = R_wc*V_bc;
    V_wc = euler_to_R(latest_mean(3),latest_mean(4), latest_mean(5))*V_bc;
    visualizeVelocity(Eigen::Vector3d(latest_mean(0), latest_mean(1), latest_mean(2)), -V_wc, 0, Eigen::Vector3d(0,1,0), opflow_vel_pub);

    if (last_tag_count != tag_count)
    {
        //FUSE TAG_DETECTOR & OPTICAL_FLOW
        MatrixXd C = MatrixXd::Identity(9,15);
        MatrixXd W = MatrixXd::Identity(9,9);
        MatrixXd K = latest_var*C.transpose()*(C*latest_var*C.transpose() + W*Rt*W.transpose()).inverse();
        state_measured(6) = -V_wc(0); 
        state_measured(7) = -V_wc(1);
        for(int i=3;i<6;i++){
            if(abs(last_rpy(i-3)-state_measured(i))>M_PI){
                if(state_measured(i)<0) state_measured(i) = 2 * M_PI + state_measured(i);
                else state_measured(i) = -2* M_PI + state_measured(i);   
            }         
        }
        last_rpy(0) = state_measured(3);
        last_rpy(1) = state_measured(4);
        last_rpy(2) = state_measured(5);
    
        latest_mean = latest_mean + K*(C*state_measured - C*latest_mean);
        latest_var = latest_var - K*C*latest_var;

        last_tag_count = tag_count;
    } 
    else
    {
        //FUSE OPTICAL_FLOW & HEIGHT ONLY
        MatrixXd C = MatrixXd::Zero(4,15); 
        C(0,2) = 1;
        C.block<3, 3>(1, 6) = I3;
        MatrixXd W = MatrixXd::Identity(4,4);
        MatrixXd K = latest_var*C.transpose()*(C*latest_var*C.transpose() + W*Rt2*W.transpose()).inverse();
        state_measured(2) = height;
        state_measured(6) = -V_wc(0); 
        state_measured(7) = -V_wc(1);
        latest_mean = latest_mean + K*(C*state_measured - C*latest_mean);
        latest_var = latest_var - K*C*latest_var;
    }
}




int main(int argc, char **argv)
{
    
    ros::init(argc, argv, "ekf");
    ros::NodeHandle n("~");
    ros::Duration(1).sleep();
    ros::Subscriber s1 = n.subscribe("imu", 1000, imu_callback);
    ros::Subscriber s2 = n.subscribe("tag_odom", 1, odom_callback_tag);
    ros::Subscriber s3 = n.subscribe("/optical_flow/velocity", 1, odom_callback_opflow);
    ros::Subscriber sub_vicon = n.subscribe("/uwb_vicon_odom", 10, viconCallback);
    ekf_odom_pub = n.advertise<nav_msgs::Odometry>("ekf_odom", 100);
    ref_odom_pub = n.advertise<nav_msgs::Odometry>("ref_odom", 100);
    opflow_vel_pub = n.advertise<visualization_msgs::Marker>("opflow_vel", 1, true);
    vicon_vel_pub = n.advertise<visualization_msgs::Marker>("vicon_vel", 1, true);
    ekf_vel_pub = n.advertise<visualization_msgs::Marker>("ekf_vel", 1, true);
    // Q imu covariance matrix; Rt visual odomtry covariance matrix
    // You should also tune these parameters
    Q.topLeftCorner(6, 6) = 0.01 * Q.topLeftCorner(6, 6);
    Q.bottomRightCorner(6, 6) = 0.01 * Q.bottomRightCorner(6, 6);
    Rt.block<3, 3>(0, 0) = 0.1 * Rt.block<3, 3>(0, 0);
    Rt.block<3, 3>(3, 3) = 0.1 * Rt.block<3, 3>(3, 3);
    Rt(5,5) = 0.1 * Rt(5,5);
    Rt.block<3, 3>(6, 6) = 0.1 * Rt.block<3, 3>(6, 6);
    Rt(8,8) = 0.1;


    Rt2(0,0) = 0.1;
    Rt2.block<3, 3>(1, 1) = Rt.block<3, 3>(6, 6);    
    ros::spin();
}
