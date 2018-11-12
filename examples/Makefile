nbExec/evol.ipynb: nb/evol.ipynb
	jupyter nbconvert  --ExecutePreprocessor.timeout=80 --to notebook --execute  $^ --output ../$@

nb/evol.ipynb: py/evol.py
	jupytext --to notebook  $^; mv  py/$(notdir $@) $@
