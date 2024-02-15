## Instructions
1. Install packages: ```./toolchains/ubuntu/install_packages.sh```
2. Configure cmake: ```cmake . -B build --preset gcc -DCMAKE_BUILD_TYPE=Release```
3. Build: ```cmake --build build --config Release```
4. Run: ```./build/src/fm_demod_tuner```