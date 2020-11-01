
# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: build
	@cd build ; make --no-print-directory -j8

build: CMakeLists.txt
	@rm -rf build
	@mkdir build
	@cd build ; cmake ..

.PHONY: debug
debug:
	@rm -rf build
	@mkdir build
	@cd build ; cmake -DDEBUG_MODE=ON ..
	@make --no-print-directory

.PHONY: release
release:
	@rm -rf build
	@mkdir build
	@cd build ; cmake ..
	@make --no-print-directory

.PHONY: clean
clean:
	@rm -rf build
	@rm -rf test

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find src/ -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================

.PHONY: test
test: debug
	@rm -rf test
	@mkdir test
	@echo "\n\e[35m\e[1m== WF run ===============================================================\e[0m"
	build/WrightFisher --sequence_size 100 --mutation_rate 0.001 --selection_coefficient 0.1 --population_size 100 --time 1000 --output test/wf
