## magicgui tests

To get current z from napari:  
```python
viewer.dims.point[0]
``` 

```python

@magicgui(auto_call=True,
          currentZ={'widget_type': 'Slider', 'min': viewer.dims.range[0][0], 'max': viewer.dims.range[0][1], 'value':viewer.dims.point[0]},
          loZ={'widget_type': 'Slider', 'min': viewer.dims.range[0][0], 'max': viewer.dims.range[0][1], 'value':viewer.dims.range[0][0]},
          hiZ={'widget_type': 'Slider', 'min': viewer.dims.range[0][0], 'max': viewer.dims.range[0][1], 'value':viewer.dims.range[0][1]},
          set_loZ = {'widget_type': 'PushButton', 'text':'set low Z'}
         )
def zlims_Widget(currentZ: int, loZ: int, hiZ: int, set_loZ: bool):
    viewer.dims.set_point(0,currentZ)
    

viewer.window.add_dock_widget(zlims_Widget)
```

#### Another attempt
> Challenging to bind values from 2 widgets. 

```python
@magicgui(auto_call=True,
          currentZ={'widget_type': 'Slider',
                    'min': viewer.dims.range[0][0], 'max': viewer.dims.range[0][1],
                    'value':viewer.dims.point[0]},
          # loZ={'widget_type': 'Slider',
          #      'min': viewer.dims.range[0][0], 'max': viewer.dims.range[0][1],
          #      'value':viewer.dims.range[0][0]},
          # hiZ={'widget_type': 'Slider',
          #      'min': viewer.dims.range[0][0], 'max': viewer.dims.range[0][1],
          #      'value':viewer.dims.range[0][1]},
          # setZ={
          #   "widget_type": "RadioButtons",
          #   "orientation": "horizontal",
          #   "choices": [("loZ", 1), ("hiZ", 2), ("setZ",3)]},
         )
def zlims_Widget(currentZ: int): #, loZ: int, hiZ: int, setZ = 1):
    viewer.dims.set_point(0,currentZ)
    # if setZ==1:
    #     loZ = currentZ
    # elif setZ ==2:
    #     hiZ = currentZ
    return currentZ
    
@magicgui(auto_call=True,
          loZ={'widget_type': 'Slider',
               'min': viewer.dims.range[0][0], 'max': viewer.dims.range[0][1],
               'value':viewer.dims.range[0][0]},
         )
def loZ_Widget(loZ: int):
    # loZ=zlims_Widget()# currentZ
    loZ=currentZ
    
@zlims_Widget.changed.connect
def _on_currentZ_changed(currentZ: int):
    # print(f"changed to {currentZ}")
    loZ_Widget(currentZ)
    
viewer.window.add_dock_widget([zlims_Widget, loZ_Widget])
```