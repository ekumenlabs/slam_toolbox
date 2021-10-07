#!/bin/bash

apt-get update
source /opt/ros/foxy/setup.bash
cd /colcon_ws
rosdep update
rosdep install -y -r --from-paths src --ignore-src --rosdistro=foxy -y