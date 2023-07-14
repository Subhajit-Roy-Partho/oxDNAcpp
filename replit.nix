{ pkgs }: {
	deps = [
		pkgs.crystal_0_35.bin
  pkgs.clang_12
		pkgs.ccls
		pkgs.gdb
		pkgs.gnumake
    pkgs.eigen
    pkgs.gsl
	];
}