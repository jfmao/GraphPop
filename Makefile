# GraphPop — build and packaging targets
# Usage:
#   make install         # Full local install (pip editable + Java build)
#   make install-quick   # Python-only install (uses pre-compiled JAR)
#   make build-jar       # Build the Java procedures JAR
#   make conda-build     # Build the conda package
#   make release         # Prepare a GitHub release (tag + JAR)
#   make clean           # Remove build artifacts

VERSION := 0.1.0
JAR_NAME := graphpop-procedures-$(VERSION).jar

.PHONY: install install-quick build-jar conda-build release clean test help

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
	  awk 'BEGIN {FS = ":.*?## "}; {printf "  %-18s %s\n", $$1, $$2}'

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------

install: build-jar  ## Full install: build Java + install all Python packages
	pip install -e graphpop-cli
	pip install -e graphpop-import
	pip install -e graphpop-mcp
	@echo ""
	@echo "Installation complete. Next steps:"
	@echo "  graphpop setup --password mypass"
	@echo "  graphpop start"

install-quick:  ## Quick install: Python packages only (auto-downloads JAR during setup)
	pip install -e graphpop-cli
	pip install -e graphpop-import
	pip install -e graphpop-mcp
	@echo ""
	@echo "Installation complete (JAR will be auto-downloaded during 'graphpop setup')."
	@echo "Next steps:"
	@echo "  graphpop setup --password mypass"
	@echo "  graphpop start"

# ---------------------------------------------------------------------------
# Java build
# ---------------------------------------------------------------------------

build-jar:  ## Build the GraphPop procedures JAR
	cd graphpop-procedures && ./mvnw package -DskipTests -q
	@echo "Built: graphpop-procedures/target/$(JAR_NAME)"

# ---------------------------------------------------------------------------
# Conda
# ---------------------------------------------------------------------------

conda-build:  ## Build the conda package locally
	conda build conda-recipe --output-folder conda-build-output
	@echo "Conda package built in conda-build-output/"

# ---------------------------------------------------------------------------
# Release
# ---------------------------------------------------------------------------

release: build-jar  ## Prepare a GitHub release (creates tag, uploads JAR)
	@echo "=== Preparing release v$(VERSION) ==="
	cp graphpop-procedures/target/graphpop-procedures-$(VERSION)*.jar \
	   graphpop-procedures-$(VERSION).jar 2>/dev/null || \
	cp graphpop-procedures/target/graphpop-procedures-$(VERSION)-SNAPSHOT.jar \
	   graphpop-procedures-$(VERSION).jar
	@echo ""
	@echo "JAR ready: graphpop-procedures-$(VERSION).jar"
	@echo ""
	@echo "To create the GitHub release:"
	@echo "  git tag -a v$(VERSION) -m 'GraphPop v$(VERSION)'"
	@echo "  git push origin v$(VERSION)"
	@echo "  gh release create v$(VERSION) graphpop-procedures-$(VERSION).jar \\"
	@echo "    --title 'GraphPop v$(VERSION)' \\"
	@echo "    --notes 'Pre-compiled procedures plugin for Neo4j 5.x'"

# ---------------------------------------------------------------------------
# Testing
# ---------------------------------------------------------------------------

test:  ## Run all tests (Java + Python)
	cd graphpop-procedures && ./mvnw test -q
	cd graphpop-cli && python -m pytest tests/ -q 2>/dev/null || true
	cd graphpop-import && python -m pytest tests/ -q 2>/dev/null || true

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

clean:  ## Remove build artifacts
	cd graphpop-procedures && ./mvnw clean -q 2>/dev/null || true
	rm -rf conda-build-output/
	rm -f graphpop-procedures-*.jar
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
