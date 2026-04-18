const std = @import("std");

const ExampleSpec = struct {
    run_args: []const []const u8,
    summary: []const u8,
    run_step: []const u8,
    run_doc: []const u8,
    physics_source: []const u8,
    physics_module: []const u8,
};

const example_specs = [_]ExampleSpec{
    .{
        .run_args = &.{"maxwell"},
        .summary = "Maxwell on 2D or 3D meshes (`--dim 2|3`)",
        .run_step = "run-maxwell",
        .run_doc = "Run the Maxwell example suite",
        .physics_source = "examples/maxwell/root.zig",
        .physics_module = "maxwell",
    },
    .{
        .run_args = &.{"euler"},
        .summary = "Incompressible Euler in 2D or 3D (`--dim 2|3`)",
        .run_step = "run-euler",
        .run_doc = "Run the Euler example suite",
        .physics_source = "examples/euler/root.zig",
        .physics_module = "euler",
    },
    .{
        .run_args = &.{ "diffusion", "--surface", "plane" },
        .summary = "Scalar diffusion on a plane or sphere",
        .run_step = "run-diffusion",
        .run_doc = "Run the diffusion example suite",
        .physics_source = "examples/diffusion/root.zig",
        .physics_module = "diffusion",
    },
};

fn createPhysicsModule(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    mod: *std.Build.Module,
    examples_common_mod: *std.Build.Module,
    source: []const u8,
) *std.Build.Module {
    return b.createModule(.{
        .root_source_file = b.path(source),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "flux", .module = mod },
            .{ .name = "examples_common", .module = examples_common_mod },
        },
    });
}

fn emitExampleTestRoot(
    b: *std.Build,
    write_files: *std.Build.Step.WriteFile,
    spec: ExampleSpec,
) std.Build.LazyPath {
    var source = std.ArrayListUnmanaged(u8){};
    defer source.deinit(b.allocator);
    const writer = source.writer(b.allocator);

    writer.writeAll("const std = @import(\"std\");\n") catch unreachable;
    writer.print("const physics = @import(\"{s}\");\n", .{spec.physics_module}) catch unreachable;
    writer.writeAll(
        \\
        \\
        \\test {
        \\    _ = physics;
        \\}
        \\
    ) catch unreachable;

    return write_files.add(
        b.fmt("generated/tests/{s}.zig", .{spec.physics_module}),
        source.items,
    );
}

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

    // -- examples_common module --
    // Shared infrastructure consumed by every example: CLI parser, snapshot
    // writer, progress display. Lives under examples/common/ and depends on
    // the flux library.
    const examples_common_mod = b.createModule(.{
        .root_source_file = b.path("examples/common/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "flux", .module = mod },
        },
    });
    const examples_common_tests = b.addTest(.{
        .root_module = examples_common_mod,
    });
    const run_examples_common_tests = b.addRunArtifact(examples_common_tests);
    const new_maxwell_mod = b.createModule(.{
        .root_source_file = b.path("examples/new_maxwell/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "flux", .module = mod },
            .{ .name = "examples_common", .module = examples_common_mod },
        },
    });
    const new_diffusion_mod = b.createModule(.{
        .root_source_file = b.path("examples/new_diffusion/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "flux", .module = mod },
            .{ .name = "examples_common", .module = examples_common_mod },
        },
    });
    const new_euler_mod = b.createModule(.{
        .root_source_file = b.path("examples/new_euler/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "flux", .module = mod },
            .{ .name = "examples_common", .module = examples_common_mod },
        },
    });
    const new_diffusion_tests = b.addTest(.{
        .root_module = new_diffusion_mod,
    });
    const run_new_diffusion_tests = b.addRunArtifact(new_diffusion_tests);
    const new_euler_tests = b.addTest(.{
        .root_module = new_euler_mod,
    });
    const run_new_euler_tests = b.addRunArtifact(new_euler_tests);
    const new_cli_commands_mod = b.createModule(.{
        .root_source_file = b.path("examples/new_cli/commands.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "new_maxwell", .module = new_maxwell_mod },
        },
    });
    const new_cli_exe = b.addExecutable(.{
        .name = "flux-new-cli",
        .root_module = b.createModule(.{
            .root_source_file = b.path("examples/new_cli/root.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "new_cli_commands", .module = new_cli_commands_mod },
            },
        }),
    });
    b.installArtifact(new_cli_exe);
    const new_cli_step = b.step("new-cli", "Build the fresh new example CLI");
    new_cli_step.dependOn(&new_cli_exe.step);
    new_cli_step.dependOn(b.getInstallStep());
    const run_new_cli_cmd = b.addRunArtifact(new_cli_exe);
    run_new_cli_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_new_cli_cmd.addArgs(args);
    }
    const run_new_cli_step = b.step("run-new-cli", "Run the fresh new example CLI");
    run_new_cli_step.dependOn(&run_new_cli_cmd.step);
    const new_cli_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("examples/new_cli/root.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "new_cli_commands", .module = new_cli_commands_mod },
            },
        }),
    });
    const run_new_cli_tests = b.addRunArtifact(new_cli_tests);
    const bench_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "flux", .module = mod },
                .{ .name = "maxwell_example", .module = b.createModule(.{
                    .root_source_file = b.path("examples/maxwell/root.zig"),
                    .target = target,
                    .optimize = optimize,
                    .imports = &.{
                        .{ .name = "flux", .module = mod },
                        .{ .name = "examples_common", .module = examples_common_mod },
                    },
                }) },
                .{ .name = "diffusion_sphere_example", .module = b.createModule(.{
                    .root_source_file = b.path("examples/diffusion/sphere.zig"),
                    .target = target,
                    .optimize = optimize,
                    .imports = &.{
                        .{ .name = "flux", .module = mod },
                        .{ .name = "examples_common", .module = examples_common_mod },
                    },
                }) },
            },
        }),
    });
    const run_bench_tests = b.addRunArtifact(bench_tests);

    const generated_sources = b.addWriteFiles();

    var physics_modules: [example_specs.len]*std.Build.Module = undefined;
    inline for (example_specs, 0..) |spec, idx| {
        physics_modules[idx] = createPhysicsModule(
            b,
            target,
            optimize,
            mod,
            examples_common_mod,
            spec.physics_source,
        );
    }

    var command_imports: [example_specs.len + 1]std.Build.Module.Import = undefined;
    command_imports[0] = .{ .name = "examples_common", .module = examples_common_mod };
    inline for (example_specs, 0..) |spec, idx| {
        command_imports[idx + 1] = .{
            .name = spec.physics_module,
            .module = physics_modules[idx],
        };
    }
    const example_commands_mod = b.createModule(.{
        .root_source_file = b.path("examples/commands.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &command_imports,
    });
    const example_app_mod = b.createModule(.{
        .root_source_file = b.path("examples/app.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "example_commands", .module = example_commands_mod },
        },
    });
    exe.root_module.addImport("example_app", example_app_mod);

    // -- flux-examples umbrella binary --
    // A single executable exposes every example as a subcommand. The build
    // graph owns the canonical example manifest for build/run/test steps,
    // while the runtime surface is the shared example app module.
    const examples_exe = b.addExecutable(.{
        .name = "flux-examples",
        .root_module = b.createModule(.{
            .root_source_file = b.path("examples/app.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "example_commands", .module = example_commands_mod },
            },
        }),
    });
    b.installArtifact(examples_exe);

    // `zig build examples` builds the umbrella binary without running it,
    // satisfying the issue acceptance criterion that all examples compile.
    const examples_step = b.step("examples", "Build the flux-examples umbrella binary");
    examples_step.dependOn(&examples_exe.step);
    examples_step.dependOn(b.getInstallStep());

    // Convenience run steps. Each `run-<name>` is a thin wrapper that runs
    // `flux-examples <name>` with the user's `--` arguments appended. The
    // build step name uses hyphens for parity with the subcommand names.
    for (example_specs) |spec| {
        const run_cmd_step = b.addRunArtifact(examples_exe);
        run_cmd_step.step.dependOn(b.getInstallStep());
        run_cmd_step.addArgs(spec.run_args);
        if (b.args) |args| {
            run_cmd_step.addArgs(args);
        }
        const top = b.step(spec.run_step, spec.run_doc);
        top.dependOn(&run_cmd_step.step);
    }
    // Each example physics module is registered once as a named module.
    // Dedicated generated test roots import those modules by name and
    // explicitly ref all decls, so adding a new `@import()` inside an example
    // can no longer silently expand its test suite.
    var example_run_test_steps: [example_specs.len + 1]*std.Build.Step.Run = undefined;
    inline for (example_specs, 0..) |spec, idx| {
        var test_imports: [1]std.Build.Module.Import = undefined;
        test_imports[0] = .{ .name = spec.physics_module, .module = physics_modules[idx] };
        const test_root = emitExampleTestRoot(b, generated_sources, spec);
        const test_mod = b.createModule(.{
            .root_source_file = test_root,
            .target = target,
            .optimize = optimize,
            .imports = &test_imports,
        });
        const test_exe = b.addTest(.{ .root_module = test_mod });
        example_run_test_steps[idx] = b.addRunArtifact(test_exe);
    }

    // The capstone target imports each physics module by name so its test
    // root only runs the five acceptance tests, not the per-example test
    // suites already covered by the modules above.
    {
        var acceptance_imports: [physics_modules.len + 2]std.Build.Module.Import = undefined;
        acceptance_imports[0] = .{ .name = "flux", .module = mod };
        acceptance_imports[1] = .{ .name = "examples_common", .module = examples_common_mod };
        inline for (example_specs, 0..) |spec, idx| {
            acceptance_imports[2 + idx] = .{ .name = spec.physics_module, .module = physics_modules[idx] };
        }
        const acceptance_mod = b.createModule(.{
            .root_source_file = b.path("examples/acceptance.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &acceptance_imports,
        });
        const acceptance_test_exe = b.addTest(.{ .root_module = acceptance_mod });
        example_run_test_steps[example_specs.len] = b.addRunArtifact(acceptance_test_exe);
    }

    // A top level step for running all tests. dependOn can be called multiple
    // times and since the two run steps do not depend on one another, this will
    // make the two of them run in parallel.
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_mod_tests.step);
    test_step.dependOn(&run_exe_tests.step);
    test_step.dependOn(&run_examples_common_tests.step);
    test_step.dependOn(&run_new_diffusion_tests.step);
    test_step.dependOn(&run_new_euler_tests.step);
    test_step.dependOn(&run_new_cli_tests.step);
    test_step.dependOn(&run_bench_tests.step);
    for (example_run_test_steps) |run_step_ptr| {
        test_step.dependOn(&run_step_ptr.step);
    }

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
                    .root_source_file = b.path("examples/maxwell/root.zig"),
                    .target = target,
                    .optimize = .ReleaseFast,
                    .imports = &.{
                        .{ .name = "flux", .module = mod },
                        .{ .name = "examples_common", .module = examples_common_mod },
                    },
                }) },
                .{ .name = "diffusion_sphere_example", .module = b.createModule(.{
                    .root_source_file = b.path("examples/diffusion/sphere.zig"),
                    .target = target,
                    .optimize = .ReleaseFast,
                    .imports = &.{
                        .{ .name = "flux", .module = mod },
                        .{ .name = "examples_common", .module = examples_common_mod },
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
    check_step.dependOn(&examples_exe.step);
    check_step.dependOn(&new_cli_exe.step);
    check_step.dependOn(&examples_common_tests.step);
    check_step.dependOn(&bench_tests.step);

    // -- ci step --
    // Runs build + test + fmt in one command: `zig build ci`
    const ci_step = b.step("ci", "Run all CI checks (build + test + fmt)");
    ci_step.dependOn(b.getInstallStep());
    ci_step.dependOn(&run_mod_tests.step);
    ci_step.dependOn(&run_exe_tests.step);
    ci_step.dependOn(&run_examples_common_tests.step);
    ci_step.dependOn(&run_new_cli_tests.step);
    ci_step.dependOn(&run_bench_tests.step);
    ci_step.dependOn(&examples_exe.step);
    ci_step.dependOn(&new_cli_exe.step);
    for (example_run_test_steps) |run_step_ptr| {
        ci_step.dependOn(&run_step_ptr.step);
    }
    ci_step.dependOn(&fmt_cmd.step);
}
