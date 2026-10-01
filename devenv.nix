{ pkgs, lib, config, inputs, ... }:

# Runtime libraries for the eframe/egui renderer.
#
# - winit talks to the Wayland compositor and falls back to X11 through x11-dl,
#   which dlopens libX11/libXcursor/libXi/libXrandr at runtime.
# - glutin (the glow backend) needs a GL/EGL implementation.
# - the `wgpu` feature of eframe needs the Vulkan loader.
#
# The Rust channel is pinned to a dated nightly. The crate no longer needs
# nightly (the last `#![feature(iter_array_chunks)]` gate went away when the
# renderer stopped using `Iterator::array_chunks`), so this can be switched to
# `channel = "stable"` with no source change once you want to.
let
  waylandLibs = with pkgs; [
    wayland
    libxkbcommon
  ];

  x11Libs = with pkgs; [
    xorg.libX11
    xorg.libXcursor
    xorg.libXi
    xorg.libXrandr
    xorg.libXrender
  ];

  # libEGL and libGLX now live in the libglvnd output set; plain `libEGL` is gone.
  glLibs = with pkgs; [
    libGL
    libglvnd
  ];

  gpuLibs = with pkgs; [
    vulkan-loader
  ];

  appLibs = waylandLibs ++ x11Libs ++ glLibs ++ gpuLibs;
in
{
  # https://devenv.sh/packages/
  packages = appLibs;

  env = {
    GREET = "devenv";
    LD_LIBRARY_PATH = "${lib.makeLibraryPath appLibs}";
    INCLUDE_DIRECTORIES = "${lib.makeIncludePath appLibs}";
  };

  # https://devenv.sh/languages/
  languages.rust = {
    enable = true;
    # rust-overlay (see devenv.yaml) supplies dated toolchains. `version` only has
    # an effect while `channel` is not "nixpkgs".
    channel = "nightly";
    version = "2026-09-27";
  };

  # https://devenv.sh/scripts/
  scripts.hello.exec = ''
    echo hello from $GREET
  '';

  # https://devenv.sh/tasks/
  # tasks = {
  #   "light_transport:check".exec = "cargo check";
  #   "devenv:enterShell".after = [ "light_transport:check" ];
  # };

  # https://devenv.sh/services/
  # services.postgres.enable = true;

  # https://devenv.sh/processes/
  # processes.dev.exec = "${lib.getExe pkgs.watchexec} -n -- cargo run";

  # https://devenv.sh/git-hooks/
  # git-hooks.hooks.shellcheck.enable = true;

  enterShell = ''
    echo "hello from $GREET"
    cargo --version
    git --version
  '';

  enterTest = ''
    echo "Running tests"
    echo "Verifying the pinned toolchain"
    cargo --version | grep --color=auto "nightly"

    echo "Type-checking the crate"
    cargo check --quiet
  '';

  # See full reference at https://devenv.sh/reference/options/
}
