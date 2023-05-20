echo -e "blacklist nouveau" | tee -a /etc/modprobe.d/blacklist.conf
echo -e "blacklist nouveau\noptions nouveau modeset=0" > /etc/modprobe.d/blacklist-nouveau.conf
mv /boot/initramfs-$(uname -r).img /boot/initramfs-$(uname -r).img.bak
update-initramfs -u
dracut -v /boot/initramfs-$(uname -r).img $(uname -r)

chmod +x NVIDIA-Linux-x86_64-465.24.02.run
./NVIDIA-Linux-x86_64-465.24.02.run --ui=none --no-questions --accept-license --disable-nouveau --no-cc-version-check --install-libglvnd
sleep 1
chmod +x cuda_11.3.0_465.19.01_linux.run
./cuda_11.3.0_465.19.01_linux.run --override --silent --toolkit
mkdir -p install && mv NVIDIA-Linux-x86_64-465.24.02.run cuda_11.3.0_465.19.01_linux.run install

cat /proc/driver/nvidia/version
lsmod | grep nvidia
rmmod nvidia_drm
rmmod nvidia_uvm
rmmod nvidia_modeset
rmmod nvidia
