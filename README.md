# MinimumSnap-Trajectory-Generation
---
Implementation of minimum snap trajectory generation for quadrotors in ROS

## 1. Prerequisites
+ Ubuntu 64-bit 16.04
+ ROS Kinetic
+ Eigen

## 2.Build on ROS
+ Clone this repository to your catkin workspace and catkin_make.
```
  cd ${YOUR_WORKSPACE_PATH}/src
  git clone https://github.com/Mesywang/MinimumSnap-Trajectory-Generation.git
  cd ../
  catkin_make
```
## 3. Run the Demo

+ Launch the package.
```
  roslaunch waypoint_trajectory_generator test.launch  
```

## 4. The result
　　Normally,you can see '3D Nav Goal' in the rviz tools panel. To set some waypoints, click the '3D Nav Goal', Then click and hold both the left and right mouse buttons to select (x,y), and move the mouse to change z. you need to set a waypoint z < 0 as the end to confirm all waypoints. After these, you will see a path marked in green and a minimum snap trajectory marked in red, as displayed below:

<div align=center>
	<img src="./img/minimum snap trajectory.gif" width = "620" height = "420" >
</div>

