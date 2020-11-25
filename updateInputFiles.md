Adapting old input files to latest version of Gkeyll
----------------------------------------------------

In this file we list a few changes necessary to run old input files with the latest version of Gkeyll. It is also helpful to look at the input files in the `gkyl/Regression` folder of the gkyl repo as examples of up-to-date input file format.

1. New way of specifying the App

See [this commit message](https://github.com/ammarhakim/gkyl/commit/7ff8debfa271959a7c73f6e1184c3837051a3ecc). Basically
older input files specified the App at the top of the file in a variety of ways, including for example

```lua
  local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell
```

Even older input files simply asked for the plasma App without specifying a model (e.g. Vlasov or Gyrokinetic), through

```lua
  local Plasma = require("App.PlasmaOnCartGrid")
```

and later narrowed the scope with calls to `Plasma.VlasovSpecies.module`, where `module` could be on of many App components (e.g. geometry, fields, etc). Now Gkeyll expects you to specify the model as in the first of these two, and add two parenthesis at the end. So please use

```lua
  local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
```
and similarly for other models (e.g. gyrokinetics, fluids).

2. Specifying components of the App

The above change (1.) also means that one should use generic calls to components of an App (fields, species, sources, etc) and not species or app-specific ones. This means that older gyrokinetic input files invoked the time varying fields through

```lua
  field = Plasma.GkField { ...
```

and the external magnetic field/geometry through

```lua
  funcField = Plasma.GkGeometry { ...
```

Provided that you follow (1.) above and specify the app as `local Plasma = require("App.PlasmaOnCartGrid").Gyrokinetic()`, one should now invoke the field and the geometry using

```lua
  field = Plasma.Field { ...

  funcField = Plasma.Geometry { ...
```

The same applies to other components (sources, species, initialization functions, etc).
