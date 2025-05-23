---
layout: post
title: Installation notes for ABINIT
date: 2019-11-29 09:20:13
description: 基于Intel2015+impi编译器编译并行版Abinit软件的方法说明
tags: compile
categories: Abinit
related_posts: false
mathjax: true
---



编译环境：Intel2015+impi

# 准备工作

* 在[Abinit](https://www.abinit.org/)官网下载最新版安装包abinit-8.10.3。
* 下载附加程序包[这里](https://www.abinit.org/fallbacks)：LibXC 3.0.0，NetCDF 4.1.1和LAPACK for Abinit ≥ 6.10。

# 安装过程

* 上传安装包与附加程序包到指定安装路径，此处以~/software/为例。
* 解压缩并在abinit-8.10.3中新建build与tarballs文件夹。
``` bash
tar zxvf abinit-8.10.3.tar.gz && cd abinit-8.10.3 && mkdir build tarballs
```
* 将LibXC 3.0.0，NetCDF 4.1.1和LAPACK for Abinit ≥ 6.10拷贝到tarballs中。
* 进入build文件夹中，并查看计算机名称。
``` bash
hostname
```
* 在build中建立hostname.ac文件，内容如下：
```
# ================================================================
# Configuration file for ABINIT 8 compilation on COBALT
# tested for Intel2015 + impi
#
# ================================================================
#
 FC="mpiifort"
 CC="mpiicc"
 CXX="mpicxx"
#
 enable_mpi="yes"
 enable_openmp="yes"
#
 with_linalg_flavor="mkl+scalapack"
 with_linalg_libs=${SCALAPACK_LDFLAGS}
#
 with_fft_flavor="fftw3"
 with_fft_incs="-I${MKL_INCDIR}"
 with_fft_libs=${MKL_LDFLAGS}
#
 with_trio_flavor="netcdf"
 with_dft_flavor="libxc"
```
* 在build中执行configure如下：
```
../configure --with-tardir=~/software/abinit-8.10.3/tarballs !根据自己机器路径修改
```
* 如果configure过程没有报错，即可编译：
``` bash
make mj4
```

# 编译测试
在build中执行：
```
cd tests && ../../tests/runtests.py fast
```
如果测试通过，则完成编译。
编译好的可执行文件为：
``` bash
~/software/abinit-8.10.3/build/src/98_main/abinit 
```

# 参考资料
* Abinit官网：https://www.abinit.org/
* Course and Hands-on session material：https://school2019.abinit.org/course-material
