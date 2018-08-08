platform='unknown'
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
   platform='linux'
elif [[ "$unamestr" == 'Darwin' ]]; then
   platform='osx'
fi

os=${OSTYPE//[0-9.]/}

if [[ "$os" == 'linux-gnu' ]]; then
	plat='-linux-64'
elif [[ "$os" == 'darwin' ]]; then
	plat='-osx-64'
else
	plat='-unknown'
fi
echo -n $plat
