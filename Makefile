# ── Compiler ────────────────────────────────────────────────
FC        = gfortran
FFLAGS    = -O2 -Wall
FFLAGS_DBG = -g -fcheck=all -fbacktrace -Wall -Wextra

# ── Binaries & stamp file ────────────────────────────────────
ORIGINAL_BIN := original/tester
CODE_BIN     := code/tester
OUTPUT_FILE := output.txt

#MATLAB_CMD := /mnt/c/Program\ Files/MATLAB/R2025b/bin/matlab.exe -nodisplay -nosplash -r "shape"
MATLAB_CMD :=  /mnt/c/Program\ Files/MATLAB/R2025b/bin/matlab.exe -batch "shape"

.PHONY: all original shape code clean

# ── Full pipeline: original → shape → code ──────────────────
all:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory original
	@$(MAKE) --no-print-directory shape
	@$(MAKE) --no-print-directory code

# ── Step 1: Build & run original/ ───────────────────────────
original:
	@echo ">>> [1/3] Compiling original..."
	cd original && $(FC) $(FFLAGS) -o tester *.f
	@echo ">>> [1/3] Running original..."
	cd original && ./tester 2>&1 | tee ../$(OUTPUT_FILE)

# ── Step 2: Run shape script ─────────────────────────────────
shape:
	@echo ">>> [2/3] Running shape (MATLAB)..."
	$(MATLAB_CMD)

# ── Step 3: Build & run code/ ────────────────────────────────
code:
	@echo ">>> [3/3] Compiling code..."
	cd code && $(FC) $(FFLAGS) -o tester *.f
	@echo ">>> [3/3] Running code..."
	cd code && ./tester 2>&1 | tee ../$(OUTPUT_FILE)

# ── Clean ────────────────────────────────────────────────────
clean:
	@rm -f $(ORIGINAL_BIN) $(CODE_BIN) $(OUTPUT_FILE)
	@rm -f original/res_form_fin_0.4_N60.dat original/essai.dat \
	        original/exact2.dat original/interpol2.dat original/nodes10.dat
	@rm -f code/res_form_fin_0.4_N60.dat code/essai.dat \
	        code/exact2.dat code/interpol2.dat \
	        code/fort.24 code/fort.20 code/fort.30 code/fort.40
	@echo "Cleaned."
