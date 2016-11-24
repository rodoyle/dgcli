# dgcli [![Build Status](https://travis-ci.org/DeskGen/dgcli.svg?branch=master)](https://travis-ci.org/DeskGen/dgcli)
Command Line Client for DeskGen Platform  

## Installation

### Prerequisites
1. Install Python 2.7
2. Create a directory and into it download both dgparse and dgcli from github
3. In each directory, run ```sudo python setup.py develop```

### Set Up .dgrc Configuration File

```
cp dgcli/config/sample.dgrc ~/.dgrc
```
Use ```chown``` to change permissions if required.

## Upgrade

### Method 1: use git
```
cd dgcli
git fetch
git reset --hard origin/{{ target_branch }}
sudo python setup.py develop
```

### Method 2: download

1. Remove existing code with ```rm -rf dgcli```

2. Go to https://github.com/deskgen/dgcli

3. Select release branch

4. Download and unpack new release

5. Install with ```sudo python setup.py develop```

## Quick Start

Run:
```dg``` To see a list of commands

To extract features, run:

```
dg extract ~/glob/to/match/one/or/more/*.files
```
