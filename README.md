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

**(ii) Enumeration of structures from specified kekulene/clarene units**

**(iii) Enumeration of structures with a given number of rings**


### II. Detailed instructions with examples

#### 1. Kekulenes

#### 2. Clarenes

#### 3. K-infinitenes

#### 4. C-infinitenes

#### 5. CC-infinitenes
