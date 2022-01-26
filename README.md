# oxford
This directory contains the code I used for building the prototype for my dissertation project of the MSc Software Engineering.

Files:

-pipeline_streamlined.R: how to run the entire analytical pipeline from start to end given a sample gene (here 'Stat3'). This version excludes the installation of all libraries required, so if you don't have these you will have to remove the '#' signs before the lines containing code. Therefore, running pipeline_streamlined.R makes things faster every time you want to run it

-pipeline.RData: the result of running the pipeline with a sample gene (here 'Stat3'). It contains all the resulting objects. This file is provided as a sample file for running shinyInterface_local.R, which is the code that produces the actual interface of the prototype with the information gathered.

-shinyInterface_local.R: takes the result of running the entire analytical pipeline (an example of which is pipeline_RData, which can be used directly) and produces a visual interface displaying all the information.

-the folder TestReactive_forDigitalOceanCloud/ contains the implementation in the cloud (Digital Ocean)

DMS
