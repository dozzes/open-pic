#! /bin/sh
OPIC_DIR=/home/dosipyan/build_place/opic/open-pic
echo "OPIC_DIR=$OPIC_DIR"
export LD_LIBRARY_PATH=$OPIC_DIR/3rdparty/luabind-0.8.1/stage:$LD_LIBRARY_PATH
export PATH=$OPIC_DIR:$PATH
