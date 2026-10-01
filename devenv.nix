{ pkgs, lib, config, inputs, ... }:

# Runtime libraries for the eframe/egui renderer.
#
# - winit talks to the Wayland compositor and falls back to X11 through x11-dl,
#   which dlopens libX11/libXcursor/libXi/libXrandr at runtime.
# - glutin (the glow backend) needs a GL/EGL implementation.
#
# The crate builds and tests on stable (verified against 1.99.0), so the
# toolchain is not pinned to a dated nightly. Bump `version` only if a
# nightly-only feature is ever needed again.
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

  # eframe's `wgpu` feature is off, so no Vulkan loader is needed.
  appLibs = waylandLibs ++ x11Libs ++ glLibs;
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
    # rust-overlay (see devenv.yaml) supplies dated toolchains when a
    # `version` is set; with `channel = "stable"` the plain stable channel
    # is used and there is nothing to pin.
    channel = "stable";
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
    cargo --version

    echo "Type-checking the crate"
    cargo check --quiet
  '';

  # See full reference at https://devenv.sh/reference/options/
}
