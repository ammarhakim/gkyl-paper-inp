Adapting old input files to latest version of Gkeyll
----------------------------------------------------

In this file we list a few changes necessary to run old input files with the latest version of Gkeyll.

1. New way of specifying the App

See [this commit message](https://github.com/ammarhakim/gkyl/commit/7ff8debfa271959a7c73f6e1184c3837051a3ecc). Basically
older input files specified the App in a variety of ways, including for example

```lua
  local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell
```

Now Gkeyll expects you to add two parenthesis at the end. So please use

```lua
  local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
```
