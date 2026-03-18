<p align="center">
  <h1 align="center">
    <code>omkit</code>
  </h1>

  <p align="center">
    <a href="https://github.com/stjude-rust-labs/omkit/actions/workflows/CI.yml" target="_blank">
      <img alt="CI: Status" src="https://github.com/stjude-rust-labs/omkit/actions/workflows/CI.yml/badge.svg" />
    </a>
    <a href="https://crates.io/crates/omkit" target="_blank">
      <img alt="crates.io version" src="https://img.shields.io/crates/v/omkit">
    </a>
    <img alt="crates.io downloads" src="https://img.shields.io/crates/d/omkit">
    <a href="https://github.com/stjude-rust-labs/omkit/blob/main/LICENSE-APACHE" target="_blank">
      <img alt="License: Apache 2.0" src="https://img.shields.io/badge/license-Apache 2.0-blue.svg" />
    </a>
    <a href="https://github.com/stjude-rust-labs/omkit/blob/main/LICENSE-MIT" target="_blank">
      <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg" />
    </a>
  </p>

  <p align="center">
    A cross-platform command line toolkit for omics-based analysis.
    <br />
    <a href="https://github.com/stjude-rust-labs/omkit/issues/new?assignees=&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    ·
    <a href="https://github.com/stjude-rust-labs/omkit/issues/new?assignees=&title=Descriptive%20Title&labels=bug">Report Bug</a>
    ·
    ⭐ Consider starring the repo! ⭐
    <br />

  </p>
</p>

> [!WARNING]
>
> **Project Status**: 🔴 Incubating 🔴
>
> This tool is incubating and is not ready for production use—neither in terms of
> correctness nor stability. We welcome anyone to try it out and give feedback,
> but beware that the interface and behavior may change out from under you at any time.

## 🎨 Features

`omkit` is a command line toolkit for omics-based analysis built on top of
[noodles](https://github.com/zaeleus/noodles). The following features are
top-level goals of the project.

- **Coming soon.** The toolkit is under active development.

## 🖥️ Development

To bootstrap a development environment, please use the following commands.

```bash
# Clone the repository
git clone git@github.com:stjude-rust-labs/omkit.git
cd omkit

# Build the tool in release mode
cargo build --release

# Run the tool
cargo run --release -- --help
```

## 🚧️ Tests

Before submitting any pull requests, please make sure the code passes the
following checks (from the root directory).

```bash
# Run the project's tests.
cargo test --all-features

# Ensure the project doesn't have any linting warnings.
cargo clippy --all-features

# Ensure the project passes `cargo fmt`.
cargo fmt --check

# Ensure the docs build.
cargo doc
```

## 🤝 Contributing

Contributions, issues and feature requests are welcome! Feel free to check the
[issues page](https://github.com/stjude-rust-labs/omkit/issues).

## 📝 License

This project is licensed as either [Apache 2.0][license-apache] or
[MIT][license-mit] at your discretion. Additionally, please see [the
disclaimer](https://github.com/stjude-rust-labs#disclaimer) that applies to all
crates and command line tools made available by St. Jude Rust Labs.

Copyright © 2026-Present [St. Jude Children's Research Hospital](https://github.com/stjude).

[license-apache]: https://github.com/stjude-rust-labs/omkit/blob/main/LICENSE-APACHE
[license-mit]: https://github.com/stjude-rust-labs/omkit/blob/main/LICENSE-MIT
