#!/bin/bash

CONTAINER=m3rsm-image
IMAGE_NAME=m3rsm-image

SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
REPO_DIR=`readlink -f ${SCRIPTS_DIR}/../`

DOCKER_MOUNT_ARGS="-v ${REPO_DIR}/:/colcon_ws/src"

xhost +
docker run --name ${IMAGE_NAME} --privileged --rm \
    ${DOCKER_MOUNT_ARGS} \
    -e USER=$USER -e USERID=$UID \
    --net=host \
    -it ${CONTAINER}
