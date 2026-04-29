# ECOSTRESS Level 2 Surface Temperature

This will be the repository for the ECOsystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) collection 3 level 2 surface temperature data product algorithm.

The ECOSTRESS collection 3 level 2 surface temperature data product is the pre-cursor to the [Surface Biology and Geology (SBG) collection 1 level 2 surface temperature data product algorithm](https://github.com/sbg-tir/SBG-TIR-L2-LSTE).

## Building

### ✓ Recommended: Using Conda (tested and verified)

This is the **easiest and most portable approach** that works on macOS, Linux, and Windows:

```bash
# Create a conda environment with all dependencies
conda create -n ecostress -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config -y

# Activate the environment
conda activate ecostress

# Build the code
cd src
make
```

The executable will be at `src/L2_PGE`.

**Quick one-liner:**
```bash
conda create -n ecostress -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config -y && conda activate ecostress && cd src && make
```

### Manual Build

If you already have all dependencies installed:

```bash
cd src
make
```

### Build Commands Reference

```bash
# Clean old build artifacts
cd src
make clean

# Regular build
make

# Debug build (with symbols, no optimization)
DEBUG=1 make
```

### Troubleshooting

**"HDF4 not found"**
- **Solution:** Use conda—HDF4 is not available in Homebrew
- Create environment: `conda create -n ecostress -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config -y`

**"libxml2 not found"** or **"expat not found"**
- **Solution:** Use conda (see above)—these packages don't always provide pkg-config files
- Conda handles this automatically

**"pkg-config: command not found"**
- macOS: `brew install pkg-config` or use conda (included in environment)
- Ubuntu: `sudo apt-get install pkg-config`
- CentOS: `sudo yum install pkgconfig`

**Build fails with "format specifies type" warnings**
- These are harmless; the executable is still created and works
- To suppress: add `-Wno-format` to `ADD_CFLAGS` in `src/Makefile`

**On macOS with M1/M2 (ARM):**
- Conda automatically installs ARM-native binaries—no special action needed

**Still having issues?**
- Clean and rebuild: `cd src && make clean && make`
- Verify environment: `conda list | grep -E "hdf|eccodes|libxml"`
- Check pkg-config: `pkg-config --list-all | grep hdf`
