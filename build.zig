const std = @import("std");

// Although this function looks imperative, it does not perform the build
// directly and instead it mutates the build graph (`b`) that will be then
// executed by an external runner. The functions in `std.Build` implement a DSL
// for defining build steps and express dependencies between them, allowing the
// build runner to parallelize the build automatically (and the cache system to
// know when a step doesn't need to be re-run).
pub fn build(b: *std.Build) void {
    // Standard target options allow the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});
    // Standard optimization options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall. Here we do not
    // set a preferred release mode, allowing the user to decide how to optimize.
    const optimize = b.standardOptimizeOption(.{});
    // It's also possible to define more custom flags to toggle optional features
    // of this build script using `b.option()`. All defined flags (including
    // target and optimize options) will be listed when running `zig build --help`
    // in this directory.

    // This creates a module, which represents a collection of source files alongside
    // some compilation options, such as optimization mode and linked system libraries.
    // Zig modules are the preferred way of making Zig code available to consumers.
    // addModule defines a module that we intend to make available for importing
    // to our consumers. We must give it a name because a Zig package can expose
    // multiple modules and consumers will need to be able to specify which
    // module they want to access.
    const mod = b.addModule("flux", .{
        // The root source file is the "entry point" of this module. Users of
        // this module will only be able to access public declarations contained
        // in this file, which means that if you have declarations that you
        // intend to expose to consumers that were defined in other files part
        // of this module, you will have to make sure to re-export them from
        // the root file.
        .root_source_file = b.path("src/root.zig"),
        // Later on we'll use this module as the root module of a test executable
        // which requires us to specify a target.
        .target = target,
    });

    // Here we define an executable. An executable needs to have a root module
    // which needs to expose a `main` function. While we could add a main function
    // to the module defined above, it's sometimes preferable to split business
    // logic and the CLI into two separate modules.
    //
    // If your goal is to create a Zig library for others to use, consider if
    // it might benefit from also exposing a CLI tool. A parser library for a
    // data serialization format could also bundle a CLI syntax checker, for example.
    //
    // If instead your goal is to create an executable, consider if users might
    // be interested in also being able to embed the core functionality of your
    // program in their own executable in order to avoid the overhead involved in
    // subprocessing your CLI tool.
    //
    // If neither case applies to you, feel free to delete the declaration you
    // don't need and to put everything under a single module.
    const exe = b.addExecutable(.{
        .name = "flux",
        .root_module = b.createModule(.{
            // b.createModule defines a new module just like b.addModule but,
            // unlike b.addModule, it does not expose the module to consumers of
            // this package, which is why in this case we don't have to give it a name.
            .root_source_file = b.path("src/main.zig"),
            // Target and optimization levels must be explicitly wired in when
            // defining an executable or library (in the root module), and you
            // can also hardcode a specific target for an executable or library
            // definition if desireable (e.g. firmware for embedded devices).
            .target = target,
            .optimize = optimize,
            // List of modules available for import in source files part of the
            // root module.
            .imports = &.{
                // Here "flux" is the name you will use in your source code to
                // import this module (e.g. `@import("flux")`). The name is
                // repeated because you are allowed to rename your imports, which
                // can be extremely useful in case of collisions (which can happen
                // importing modules from different packages).
                .{ .name = "flux", .module = mod },
            },
        }),
    });

    // This declares intent for the executable to be installed into the
    // install prefix when running `zig build` (i.e. when executing the default
    // step). By default the install prefix is `zig-out/` but can be overridden
    // by passing `--prefix` or `-p`.
    b.installArtifact(exe);

    // This creates a top level step. Top level steps have a name and can be
    // invoked by name when running `zig build` (e.g. `zig build run`).
    // This will evaluate the `run` step rather than the default step.
    // For a top level step to actually do something, it must depend on other
    // steps (e.g. a Run step, as we will see in a moment).
    const run_step = b.step("run", "Run the app");

    // This creates a RunArtifact step in the build graph. A RunArtifact step
    // invokes an executable compiled by Zig. Steps will only be executed by the
    // runner if invoked directly by the user (in the case of top level steps)
    // or if another step depends on it, so it's up to you to define when and
    // how this Run step will be executed. In our case we want to run it when
    // the user runs `zig build run`, so we create a dependency link.
    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);

    // By making the run step depend on the default step, it will be run from the
    // installation directory rather than directly from within the cache directory.
    run_cmd.step.dependOn(b.getInstallStep());

    // This allows the user to pass arguments to the application in the build
    // command itself, like this: `zig build run -- arg1 arg2 etc`
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // Creates an executable that will run `test` blocks from the provided module.
    // Here `mod` needs to define a target, which is why earlier we made sure to
    // set the releative field.
    const mod_tests = b.addTest(.{
        .root_module = mod,
    });

    // A run step that will run the test executable.
    const run_mod_tests = b.addRunArtifact(mod_tests);

    // Creates an executable that will run `test` blocks from the executable's
    // root module. Note that test executables only test one module at a time,
    // hence why we have to create two separate ones.
    const exe_tests = b.addTest(.{
        .root_module = exe.root_module,
    });

    // A run step that will run the second test executable.
    const run_exe_tests = b.addRunArtifact(exe_tests);

    // -- example-maxwell2d step --
    // Builds and runs the 2D Maxwell electromagnetic example.
    // This example consumes flux as a library dependency, proving the
    // package API works end-to-end.
    const maxwell2d_exe = b.addExecutable(.{
        .name = "maxwell_2d",
        .root_module = b.createModule(.{
            .root_source_file = b.path("examples/maxwell_2d/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "flux", .module = mod },
            },
        }),
    });
    b.installArtifact(maxwell2d_exe);
    const maxwell2d_run = b.addRunArtifact(maxwell2d_exe);
    maxwell2d_run.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        maxwell2d_run.addArgs(args);
    }
    const maxwell2d_step = b.step("example-maxwell2d", "Run the 2D Maxwell electromagnetic example");
    maxwell2d_step.dependOn(&maxwell2d_run.step);

    // Example tests — these are integration tests that verify the physics
    // module works correctly through the public flux API.
    const maxwell2d_tests = b.addTest(.{
        .root_module = maxwell2d_exe.root_module,
    });
    const run_maxwell2d_tests = b.addRunArtifact(maxwell2d_tests);

    // -- example-maxwell3d step --
    // Builds and runs the 3D Maxwell cavity example on tetrahedral meshes.
    const maxwell3d_exe = b.addExecutable(.{
        .name = "maxwell_3d",
        .root_module = b.createModule(.{
            .root_source_file = b.path("examples/maxwell_3d/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "flux", .module = mod },
            },
        }),
    });
    b.installArtifact(maxwell3d_exe);
    const maxwell3d_run = b.addRunArtifact(maxwell3d_exe);
    maxwell3d_run.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        maxwell3d_run.addArgs(args);
    }
    const maxwell3d_step = b.step("example-maxwell3d", "Run the 3D Maxwell cavity example");
    maxwell3d_step.dependOn(&maxwell3d_run.step);

    const maxwell3d_tests = b.addTest(.{
        .root_module = maxwell3d_exe.root_module,
    });
    const run_maxwell3d_tests = b.addRunArtifact(maxwell3d_tests);

    // -- example-euler2d step --
    // Builds and runs the 2D incompressible Euler vorticity-stream example.
    const euler2d_exe = b.addExecutable(.{
        .name = "euler_2d",
        .root_module = b.createModule(.{
            .root_source_file = b.path("examples/euler_2d/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "flux", .module = mod },
            },
        }),
    });
    b.installArtifact(euler2d_exe);
    const euler2d_run = b.addRunArtifact(euler2d_exe);
    euler2d_run.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        euler2d_run.addArgs(args);
    }
    const euler2d_step = b.step("example-euler2d", "Run the 2D incompressible Euler example");
    euler2d_step.dependOn(&euler2d_run.step);

    const euler2d_tests = b.addTest(.{
        .root_module = euler2d_exe.root_module,
    });
    const run_euler2d_tests = b.addRunArtifact(euler2d_tests);

    // A top level step for running all tests. dependOn can be called multiple
    // times and since the two run steps do not depend on one another, this will
    // make the two of them run in parallel.
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_mod_tests.step);
    test_step.dependOn(&run_exe_tests.step);
    test_step.dependOn(&run_maxwell2d_tests.step);
    test_step.dependOn(&run_maxwell3d_tests.step);
    test_step.dependOn(&run_euler2d_tests.step);

    // -- docs step --
    // Generates HTML documentation from doc comments: `zig build docs`
    // Output goes to zig-out/docs/ as a static HTML site.
    const lib = b.addLibrary(.{
        .name = "flux",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/root.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const install_docs = b.addInstallDirectory(.{
        .source_dir = lib.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    });
    const docs_step = b.step("docs", "Generate API documentation");
    docs_step.dependOn(&install_docs.step);

    // -- serve-docs step --
    // Builds docs, then launches a local HTTP server and opens the browser.
    // Zig autodoc is Wasm-based and requires HTTP (browsers block Wasm from file://).
    const serve_docs_exe = b.addExecutable(.{
        .name = "serve-docs",
        .root_module = b.createModule(.{
            .root_source_file = b.path("tools/serve-docs.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const serve_docs_run = b.addRunArtifact(serve_docs_exe);
    serve_docs_run.step.dependOn(&install_docs.step);
    const serve_docs_step = b.step("serve-docs", "Build docs and serve locally at http://127.0.0.1:8080");
    serve_docs_step.dependOn(&serve_docs_run.step);

    // -- bench step --
    // Builds and runs the benchmark suite in ReleaseFast mode.
    // Pass --check to compare against committed baselines and fail on regression.
    // Pass --update to overwrite baselines.json with current results.
    const bench_exe = b.addExecutable(.{
        .name = "bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/main.zig"),
            .target = target,
            .optimize = .ReleaseFast,
            .imports = &.{
                .{ .name = "flux", .module = mod },
                .{ .name = "maxwell_example", .module = b.createModule(.{
                    .root_source_file = b.path("examples/maxwell_2d/maxwell.zig"),
                    .target = target,
                    .optimize = .ReleaseFast,
                    .imports = &.{
                        .{ .name = "flux", .module = mod },
                    },
                }) },
            },
        }),
    });
    const bench_run = b.addRunArtifact(bench_exe);
    if (b.args) |args| {
        bench_run.addArgs(args);
    }
    const bench_step = b.step("bench", "Run operator benchmarks (--check to compare baselines, --update to save)");
    bench_step.dependOn(&bench_run.step);

    // -- fmt step --
    // Runs zig fmt --check on all source files. Fails if anything is unformatted.
    const fmt_step = b.step("fmt", "Check source formatting");
    const fmt_cmd = b.addFmt(.{
        .paths = &.{ "src", "bench", "examples", "build.zig" },
        .check = true,
    });
    fmt_step.dependOn(&fmt_cmd.step);

    // -- check step --
    // Compiles the main artifacts without installing or running them.
    // This is the fast path editors can use for near-live diagnostics.
    const check_step = b.step("check", "Compile main artifacts without running");
    check_step.dependOn(&exe.step);
    check_step.dependOn(&mod_tests.step);
    check_step.dependOn(&exe_tests.step);
    check_step.dependOn(&maxwell2d_exe.step);
    check_step.dependOn(&maxwell2d_tests.step);
    check_step.dependOn(&maxwell3d_exe.step);
    check_step.dependOn(&maxwell3d_tests.step);

    // -- ci step --
    // Runs build + test + fmt in one command: `zig build ci`
    const ci_step = b.step("ci", "Run all CI checks (build + test + fmt)");
    ci_step.dependOn(b.getInstallStep());
    ci_step.dependOn(&run_mod_tests.step);
    ci_step.dependOn(&run_exe_tests.step);
    ci_step.dependOn(&run_maxwell2d_tests.step);
    ci_step.dependOn(&run_maxwell3d_tests.step);
    ci_step.dependOn(&fmt_cmd.step);
}
