#!/bin/bash

source /opt/ros/foxy/setup.bash
cd /colcon_ws/
colcon build  --cmake-args=-DCMAKE_BUILD_TYPE=Release