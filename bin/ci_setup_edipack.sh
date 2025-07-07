#Building EDIpack
#Errors
set -e
set -v

cd edipack
pwd 

mkdir -p build
cd build
pwd

cmake ..

make -j

make install


CONF_PATH=$(grep CONFIG_FILE_PATH CMakeCache.txt |cut -d "=" -f2)

echo "source $CONF_PATH/edipack_config_user.sh" >> ~/.edipack_config_user

echo -e "\e[32m EDIpack compiled and sourced \e[0m"
