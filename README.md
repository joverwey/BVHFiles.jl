# BVH

[![Build Status](https://travis-ci.com/CarlBittendorf/BVH.jl.svg?branch=master)](https://travis-ci.com/CarlBittendorf/BVH.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/CarlBittendorf/BVH.jl?svg=true)](https://ci.appveyor.com/project/CarlBittendorf/BVH-jl)
[![Coverage](https://codecov.io/gh/CarlBittendorf/BVH.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/CarlBittendorf/BVH.jl)
[![Coverage](https://coveralls.io/repos/github/CarlBittendorf/BVH.jl/badge.svg?branch=master)](https://coveralls.io/github/CarlBittendorf/BVH.jl?branch=master)

BVH.jl is a package for working with BioVisionHierarchy files, a file format which stores motion capture 
data in a hierarchical structure. 

BVH.jl uses [graphs](https://github.com/JuliaGraphs/LightGraphs.jl) to represent and manipulate the data.

## Examples

![Original small](https://user-images.githubusercontent.com/85636219/121529929-2d5e5d80-c9fd-11eb-9b88-8622ed137318.png)
![Icon small](https://user-images.githubusercontent.com/85636219/121533093-3ac91700-ca00-11eb-9d48-8eff6ab156bc.png)
![Transformed small](https://user-images.githubusercontent.com/85636219/121529954-351e0200-c9fd-11eb-8820-7458841fd94f.png)
![Icon small](https://user-images.githubusercontent.com/85636219/121533093-3ac91700-ca00-11eb-9d48-8eff6ab156bc.png)
![Result small](https://user-images.githubusercontent.com/85636219/121530043-4a932c00-c9fd-11eb-8dae-dd53e47e847d.png)



https://user-images.githubusercontent.com/85636219/121430877-f6deff00-c978-11eb-944a-ec36bfbfc072.mp4


## Features

- Load BVH files
- Calculate global positions
- Removal of joints 
- Adding joints
- Adding frames
- Loss functions
- Optimizing offsets
- Changing rotation sequences
- Scaling
- Transfer movements to different hierarchies
- Save BVH files
