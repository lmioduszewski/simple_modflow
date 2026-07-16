from . import layout_template as template

### DEFINE LAYOUT PRESETS ###

class LayoutPresets:

    default = {}
    reports = {
        'xaxis': {
            'gridcolor': 'black',
            'griddash': 'solid',
            'ticks': 'outside',
        },
        'yaxis': {
            'gridcolor': 'black',
            'griddash': 'solid',
            'ticks': 'outside',
        }
    }

    def __init__(self):

        pass

    @classmethod
    def get(cls, preset: str = 'default') -> dict:
        """
        Gets a specific preset configuration for figure Layout.

        :param preset: Name of the preset to be applied. Defaults to 'default'.
        :type preset: str
        """
        preset = preset.lower()
        valid_presets = cls.list_presets()
        if not hasattr(cls, preset):
            raise ValueError(f"Invalid preset name: {preset}. \n"
                             f"Choose from: {valid_presets}")

        return getattr(cls, preset)

    @classmethod
    def list_presets(cls) -> list[str]:
        """
        Returns a list of valid preset names.
        """
        return [
            key for key in dir(cls)
            if not key.startswith("_") and isinstance(getattr(cls, key), dict)
        ]

    @classmethod
    def set(cls, preset: str):
        """updates the default layout with a preset and returns it"""
        preset_dict = cls.get(preset)
        updated = template.layout.update(preset_dict)
        return updated


if __name__ == "__main__":

    print(LayoutPresets.set('reports'))
