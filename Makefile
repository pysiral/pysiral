init:
	pip install -r requirements.txt
	python setup.py build_ext --inplace

test:
	python -m unittest discover -s tests -t tests