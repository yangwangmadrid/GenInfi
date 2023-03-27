# GenInfi
Generating infinitenes, kekulenes and clarenes

<img src="https://repository-images.githubusercontent.com/619128989/08bf4173-9a9c-48e2-89e0-eb4848209aeb" 
     alt="Logo of the GenInfi program" title="GenInfi" width=360 />
     
This is a handy program for generating structures of generalized infinitenes, generalized kekulenes, as well as clarenes. All constructed molecular structures are exported in Cartesian coordinates in \*.xyz files.


## Copyright and License
The Author of the Geninfi software is Yang Wang 
(yangwang@yzu.edu.cn; [orcid.org/0000-0003-2540-2199](https://orcid.org/0000-0003-2540-2199)). The Geninfi program is 
released under GNU General Public License v3 (GPLv3).

<img src="https://www.gnu.org/graphics/gplv3-or-later.png" 
     alt="Logo of GPLv3 or later" title="GPLv3 or later" />
 
 
## Disclaimer
The Geninfi software is provided as it is, with no warranties. The Author shall 
not be liable for any use derived from it. Feedbacks and bug reports are always 
welcome (yangwang@yzu.edu.cn). However, it is kindly reminded that the Author 
does not take on the responsibility of providing technical support.  
 
 
## How to Install
 
### Requirements
- python >= 3.6
- numpy >= 1.18.0


### Installation
This is a typical python program and requires no installation as long as Python 3 has been preinstalled in your computer. Upon downloading the Geninfi package, simply decompress it to any location as you like. Hereafter, we refer to the folder where the main program file, *geninfi.py*, is located as the source directory of Geninfi (denoted by `SRC_GENINFI`).


## How to Use

It is very simple to use GenInfi to generate structures of infinitenes, kekulenes, and clarenes. 

For a quick start, you can go to the `examples/` folder and run the `script*.sh` files under each subfolder for tests.

You can also refer to the detailed instructions with examples in the next section.


### I. General description of usage

#### 1. Command line arguments

The program is executed at command line, as follows:

```
python SRC_GENINFI/geninfi.py <kinf|cinf|ccinf|kek|clr> <parameters>
```
where `parameters` can be either of the following options, depending on the types of constructions:
```
- Nring  (Number of rings for the molecules to enumerate)
- h1 h2 h3 h4 h5 h6  (Six side lengths for specifying kekulenes/clarenes, or the clarene unit of CC-infinitenes)
- h1 h2 h3 h4 h5 h6 d  (Six side lengths and d-shift for specifying CC-infinitene)
- h1 h2 h3 h4 h5 h6 k1 k2 k3 k4 k5 k6  (Twelve side lengths for enumeration of K- and C-infinitenes)
- h1 h2 h3 h4 h5 h6 k1 k2 k3 k4 k5 k6 d (Twelve side lengths and d-shift for specifying K- and C-infinitenes)
```

#### 2. Types of structural constructions

**(i) Direct generation of a specified structure**

- By specifying 6 side lengths, `h1 h2 h3 h4 h5 h6`, one can construct the corresponding kekulene \[h1,h2,h3,h4,h5,h6\] or clarene <h1,h2,h3,h4,h5,h6>.

- For making a CC-infinitene <coronene|h1,h2,h3,h4,h5,h6>(d), one needs to provide an additional integer for the shift `d`.

- Likewise, given the 12 side lengths, `h1 h2 h3 h4 h5 h6 k1 k2 k3 k4 k5 k6`, and the shift `d`, the corresponding K-infinitene \[h1,h2,h3,h4,h5,h6|k1,k2,k3,k4,k5,k6\](d) or C-infinitene <h1,h2,h3,h4,h5,h6|k1,k2,k3,k4,k5,k6>(d) is obtained.

**NOTE:** All provided structural parameters must be valid for the desired molecule to construct. Otherwise, an error message will be prompted by GenInfi. 


**(ii) Enumeration of structures from specified kekulene/clarene units**

One can enumerate all possible structures of K-infinitenes/C-infinitenes that are composed of two given kekulene/clarene units by specifying the side lengths for the corresponding kekulene/clarenes. Note that a given pair of kekulene (or clarene) units can have many different combinations to form K-infinitene (or C-infinitene) structures, depending on their contacting sides, relative shifts and relative orientation.

For CC-infinitenes, similar enumeration can be performed, but one only needs to provide the side lengths of the constituting clarene unit.


**(iii) Enumeration of structures with a given number of rings**

The simplest (but most exhaustive) way to generate these macrocyclic compounds is to tell GenInfi how many benzene rings you want the molecule to have. Given the number of rings, `Nring`, the program will conduct a full structural enumeration.


#### 3. Outputs

- All detailed information is printed out to the screen during the execution of the program. This could be quite lengthy especially for a full enumeration of structures. It is thus recommended to redirect the screen output to a log file (using `python SRC_GENINFI/geninfi.py ... > log_file`).
- All generated structures are saved as Cartesian coordinates of atoms in \*.xyz files, which can be visualized by many softwares such as [JMol](https://jmol.sourceforge.net/).
- For each structure, the simple Hückel π energy is computed and given in the title line of the \*.xyz file.


### II. Detailed instructions with examples

The following examples can be found in the `examples/` folder. The execution command in each case is given in the script\*.sh file. The generated the \*.xyz output files are in the same folder.

#### 1. Kekulenes

- Generate kekulene \[1,3,2,3,1,4\]:
```
python ../../geninfi.py kek 1 3 2 3 1 4
```
**NOTE:** All side lengths of a kekulene must form a valid equiangular hexagon.

The output file `kek_R14-1_3_2_3_1_4.xyz` will be generated.

<img src="https://github.com/yangwangmadrid/GenInfi/blob/main/images/kek_R14-1_3_2_3_1_4.png" 
     alt="kekulene [1,3,2,3,1,4]" title="kekulene [1,3,2,3,1,4]" width=150 />

- Enumerate all possible \[15\]kekulene structures"
```
python ../../geninfi.py kek 15
```
As we will see, four structures are generated: `kek_R15-1_2_4_2_1_5.xyz`, `kek_R15-1_3_3_2_2_4.xyz`, `kek_R15-1_4_1_4_1_4.xyz`, and `kek_R15-2_3_2_3_2_3.xyz`.


#### 2. Clarenes

- Generate clarene <2,4,4,4,2,6>:
```
python ../../geninfi.py clr 2 4 4 4 2 6
```
**NOTE:** All side lengths of a clarene must be even numbers and form a valid equiangular hexagon.

The output file `clr_R22-2_4_4_4_2_6.xyz` will be generated.

<img src="https://github.com/yangwangmadrid/GenInfi/blob/main/images/clr_R22-2_4_4_4_2_6.png" 
     alt="clarene <2,4,4,4,2,6>" title="clarene <2,4,4,4,2,6>" width=180 />
     
- Enumerate all possible \[26\]clarene structures"
```
python ../../geninfi.py clr 26
```
Two possible isomers, `clr_R26-2_4_6_4_2_8.xyz` and `clr_R26-2_6_4_4_4_6.xyz` will be created. The structure of the first one is shown in the figure below.

<img src="https://github.com/yangwangmadrid/GenInfi/blob/main/images/clr_R26-2_4_6_4_2_8.png" 
     alt="clarene <2,4,6,4,2,8>" title="clarene <2,4,6,4,2,8>" width=180 />


#### 3. K-infinitenes

- Generate K-infinitene \[2,4,3,4,2,5|3,3,4,5,1,6\](1):
```
python ../../geninfi.py kinf 2 4 3 4 2 5 3 3 4 5 1 6 1
```
The constructed molecule contains 42 rings, as written in `kinf_R42-2_4_3_4_2_5-3_3_4_5_1_6_d1.xyz` and depicted as the following picture.

<img src="https://github.com/yangwangmadrid/GenInfi/blob/main/images/kinf_R42-2_4_3_4_2_5-3_3_4_5_1_6_d1.png" 
     alt="K-infinitene [2,4,3,4,2,5|3,3,4,5,1,6](1)" title="K-infinitene [2,4,3,4,2,5|3,3,4,5,1,6](1)" width=250 />

- Enumerate all possible K-infinitenes composed of kekulenes \[1,1,2,1,1,2\] and \[1,2,3,1,2,3\]"
```
python ../../geninfi.py kinf 1 1 2 1 1 2 1 2 3 1 2 3
```
As GenInfi outputs, there are as many as 30 enumerated structures: `kinf_R20-1_1_2_1_1_2-1_2_3_1_2_3_d-1.xyz`, ..., `kinf_R20-1_2_1_1_2_1-3_2_1_3_2_1_d0.xyz`.

- Enumerate all possible K-infinitenes containing 14 rings
```
python ../../geninfi.py kinf 14
```
Four possible isomers will be found: `kinf_R14-1_1_1_1_1_1-1_1_2_1_1_2_d-1.xyz`, `kinf_R14-1_1_1_1_1_1-1_1_2_1_1_2_d0.xyz`, `kinf_R14-1_1_1_1_1_1-1_1_2_1_1_2_d1.xyz`, and `kinf_R14-1_1_1_1_1_1-2_1_1_2_1_1_d-1.xyz`.


#### 4. C-infinitenes

- Generate C-infinitene <2,4,2,4,2,4|6,4,2,6,4,2>(-2):
```
python ../../geninfi.py cinf 2 4 2 4 2 4 6 4 2 6 4 2 -2
```
The obtained 42-ring infinitene, `cinf_R42-2_4_2_4_2_4-6_4_2_6_4_2_d-2.xyz`, shows a structure as follows.

<img src="https://github.com/yangwangmadrid/GenInfi/blob/main/images/cinf_R42-2_4_2_4_2_4-6_4_2_6_4_2_d-2.png" 
     alt="C-infinitene <2,4,2,4,2,4|6,4,2,6,4,2>(-2)" title="C-infinitene <2,4,2,4,2,4|6,4,2,6,4,2>(-2)" width=250 />

#### 5. CC-infinitenes
