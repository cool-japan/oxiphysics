# Security Policy

## Supported Versions

OxiPhysics is in its `0.1.x` series. Security fixes are applied to the **latest
published `0.1.x` patch release**; please upgrade to the most recent patch
before reporting an issue.

| Version | Supported          |
| ------- | ------------------ |
| 0.1.x (latest patch) | :white_check_mark: |
| older 0.1.x patches  | :x:                |

## Reporting a Vulnerability

**Please do not open a public GitHub issue for security problems.** Public
issues disclose the vulnerability before a fix is available.

Instead, report it privately through **GitHub's private vulnerability
reporting** on the [`cool-japan/oxiphysics`](https://github.com/cool-japan/oxiphysics)
repository:

1. Go to the repository's **Security** tab.
2. Click **Report a vulnerability**.
3. Provide as much detail as you can — affected crate/version, a description of
   the issue, and ideally a minimal reproduction or proof of concept.

This routes the report directly and privately to the maintainers.

### Scope note

OxiPhysics is a **simulation library**, not a network-facing service. Its use of
`unsafe` is deliberately confined to the **collision-dispatch** and
**GPU-backend** layers, where each `unsafe` site carries a documented soundness
contract. Reports about memory-safety issues in those layers (or anywhere else
in the workspace) are especially welcome.

### Response expectations

We aim, on a **best-effort** basis, to **acknowledge a report within a few
business days** and to keep you informed as we investigate and prepare a fix.
Timelines may vary with severity and maintainer availability. We appreciate
responsible disclosure and will credit reporters who wish to be acknowledged.
