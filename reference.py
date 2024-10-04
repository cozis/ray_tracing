"""
# Import the library using the alias "mi"
import mitsuba as mi
# Set the variant of the renderer
mi.set_variant('scalar_rgb')
# Load a scene
scene = mi.load_dict(mi.cornell_box())
# Render the scene
img = mi.render(scene)
# Write the rendered image to an EXR file
mi.Bitmap(img).write('cbox.exr')
"""
import mitsuba as mi
import cv2
import numpy as np

mi.set_variant('scalar_rgb')

# Define a scene with a sphere and a cube
scene_dict = {
    'type': 'scene',
    'integrator': {
        'type': 'path'
    },
    'sensor': {
        'type': 'perspective',
        'film': {
            'type': 'hdrfilm',
            'width': 800,
            'height': 600,
        },
        'sampler': {
            'type': 'independent',
            'sample_count': 64
        },
        'to_world': mi.ScalarTransform4f.look_at(
            origin=[10, 10, 10],     # Camera position
            target=[0, 0, 0],     # Camera looks at the origin
            up=[0, 1, 0]          # Up direction for the camera
        )
    },
    'light': {
        'type': 'point',
        'position': [0, 5, 0],
        'intensity': {'type': 'spectrum', 'value': 30.0}
    },
    'sphere': {
        'type': 'sphere',
        'center': [0, 0, 0],
        'radius': 1,
        'bsdf': {
            'type': 'diffuse',
            'reflectance': {'type': 'rgb', 'value': [0.8, 0.3, 0.3]}
        }
    },
    'cube': {
        'type': 'cube',
        'to_world': mi.ScalarTransform4f.translate([2, 0, 0]).scale(1),
        'bsdf': {
            'type': 'diffuse',
            'reflectance': {'type': 'rgb', 'value': [0.3, 0.8, 0.3]}
        }
    },
	'floor': {
        'type': 'cube',
        'to_world': mi.ScalarTransform4f.translate([0, -0.5, 0]).scale([5, 0.01, 5]),
        'bsdf': {
            'type': 'diffuse',
            'reflectance': {'type': 'rgb', 'value': [0.8, 0.8, 0.8]}  # Light grey floor
        }
    }
}

# Load and render the scene
scene = mi.load_dict(scene_dict)
img = mi.render(scene)

# Convert the rendered image to a numpy array
bitmap = mi.Bitmap(img)
image_np = np.array(bitmap.convert(mi.Bitmap.PixelFormat.RGB, mi.Struct.Type.UInt8))

# Display the image using OpenCV
cv2.imshow('Rendered Image', image_np)
cv2.waitKey(0)  # Press any key to close the window
cv2.destroyAllWindows()