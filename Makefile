.PHONY: debug release test documentation documentation-archive clean clean-docs

SWIFT_TARGET := SwiftyBoost
DOCS_DIR := Docs
DOC_ARCHIVE := $(DOCS_DIR)/$(SWIFT_TARGET).doccarchive
SYMBOL_GRAPH_DIR := .build/symbol-graphs
HOSTING_BASE_PATH ?= swiftyboost

format:
	swift format -i -r --configuration swift-format.json Sources/

debug:
	swift build

release:
	swift build -c release

test:
	swift test

documentation:
	# Build static site using the DocC SPM plugin
	swift package \
		--allow-writing-to-directory $(DOCS_DIR) \
		generate-documentation \
		--target $(SWIFT_TARGET) \
		--output-path $(DOCS_DIR) \
		--transform-for-static-hosting \
		--hosting-base-path $(HOSTING_BASE_PATH)

documentation-archive:
	# Build a .doccarchive using the DocC SPM plugin
	swift package \
		--allow-writing-to-directory $(DOCS_DIR) \
		generate-documentation \
		--target $(SWIFT_TARGET) \
		--output-path $(DOC_ARCHIVE)

clean-docs:
	rm -rf $(DOC_ARCHIVE) $(SYMBOL_GRAPH_DIR)

clean:
	swift package clean
	rm -rf $(DOCS_DIR) $(SYMBOL_GRAPH_DIR)
