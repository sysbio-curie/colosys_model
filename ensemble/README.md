
## Use with conda

First, build the conda environment : 
```
conda create -n bonesis -c potassco python=3.6 clingo jupyter notebook pandas matplotlib pyparsing pydot ginsim-python
```

Then activate it : 
```
conda activate bonesis
```

Install bonesis : 
```
pip install git+https://github.com/bioasp/bonesis.git
```


Finally, run the notebook: 
```
jupyter notebook
```


 