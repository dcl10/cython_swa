# cython_swa
Cython implementation of Smith Waterman Algorithm

## Installation
If you don't have a C/C++ compiler installed, find the right one for your OS:
- Linux: run `sudo apt-get install build-essential` in the terminal
- MacOS: get XCode from [https://developer.apple.com/](https://developer.apple.com/)
- Windows: find a suitable compiler from [https://wiki.python.org/moin/WindowsCompilers](https://wiki.python.org/moin/WindowsCompilers)

Once you have install a C/C++ compiler, run the following in the terminal/command prompt:
```shell
python -m venv venv
# Linux/MacOS
source .venv/bin/activate
# Windows
venv\Scripts\activate

python -m pip install -r requirements.txt
python setup.py build_ext --inplace

```

## Run the demo
```shell
python main.py

# Output
Python SWA took 35.515625s
Cython SWA took 1.031250s
Cython speed up of 34.44x
```