# Troubleshooting

## `apptainer` user-namespace error

```
root filesystem extraction failed: extract command failed: ERROR: Failed to create user namespace: user namespace disabled
```

Raised by `apptainer` (used inside the exitron-splicing module) on systems where the kernel has user namespaces disabled. If `exitronsplicing.activate: false` in your config, `apptainer` isn't required and this error goes away. Otherwise, ask your sysadmin to enable user namespaces, or run on a host that has them enabled.
