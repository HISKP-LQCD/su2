# yaml arguments parsing

Examples and tests of .yaml input files parsing.
This folder also contains examples of input files for the various programs:

* hmc-u1.input.example
* ... (to be completed)

## Instructions

1. ``` mkdir build ```
2. ``` cd build ```
3. Create the following 2 files:
   1. ```config.yaml```:
        ``` yaml
        p1: 4

        n1:
          x : 3
          y : 4

        p2: ciao
        ```
    2. ```build_project.sh```
        ``` bash
        cmake -D YAML_SRC=/opt/yaml-cpp/ ..; # change /opt/yaml-cpp with you path
        cmake -j 8 --build . 
        ```
4. ``` bash build_project.sh```
5. ```./yce```