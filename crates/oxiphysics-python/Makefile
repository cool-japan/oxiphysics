.PHONY: dev test build clean

dev:
	maturin develop --release

test: dev
	python -m pytest python/tests/ -v

build:
	maturin build --release --out dist/

clean:
	rm -rf dist/ target/wheels/
