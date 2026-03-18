# Release

- [ ] Update version in `omkit/Cargo.toml`.
- [ ] Update `omkit/CHANGELOG.md` with version and publication date.

  To get the changes since the last release:

  ```bash
  git log omkit-v0.1.0..HEAD --oneline -- omkit
  ```

- [ ] Run tests: `cargo test --all-features`.
- [ ] Run linting: `cargo clippy --all-features`.
- [ ] Run fmt: `cargo fmt --check`.
- [ ] Run doc: `cargo doc`.
- [ ] Stage changes: `git add omkit/Cargo.toml omkit/CHANGELOG.md`.
- [ ] Create git commit:
  ```
  git commit -m "release: bumps `omkit` version to v0.1.0"
  ```
- [ ] Create git tag:
  ```
  git tag omkit-v0.1.0
  ```
- [ ] Push release: `git push && git push --tags`.
- [ ] Publish the crate: `cargo publish --all-features -p omkit`.
- [ ] Go to the Releases page in GitHub, create a Release for this tag, and
      copy the notes from `omkit/CHANGELOG.md`.
