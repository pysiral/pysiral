init:
	pip install -e .[dev]
	python setup.py build_ext --inplace

test:
	python -m unittest discover -s tests -t tests