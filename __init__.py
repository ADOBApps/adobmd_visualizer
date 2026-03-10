"""
Author: Acxel Orozco
Date: 2026-03-09
Description: ADOBMD Visualizer - Visualize and export ADOBMD data files

This file is part of Quantum Analysis Helper.
Quantum Analysis Helper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Copyright (C) [2026] Acxel David Orozco Baldomero
"""

from pathlib import Path
from typing import Optional
import logging

from PySide6.QtWidgets import QWidget
from PySide6.QtGui import QIcon
from PySide6.QtCore import QObject

from plugins.libs.plugin_manager import PluginInfo

class ADOBMDVisualizer(QObject):
    """Main plugin class for ADOBMD Visualizer"""
    
    def __init__(self, plugin_info: PluginInfo):
        super().__init__()
        self.plugin_info = plugin_info
        self.ui_widget = None
        self.icon = None
        
        # Load icon if available
        plugin_dir = Path(__file__).parent
        icon_path = plugin_dir / "icon.png"
        if icon_path.exists():
            self.icon = QIcon(str(icon_path))
    
    def initialize(self) -> bool:
        """Initialize plugin resources"""
        try:
            logging.info(f"Initializing plugin {self.plugin_info.name}")
            return True
        except Exception as e:
            logging.error(f"Failed to initialize plugin: {e}")
            return False
    
    def get_widget(self) -> Optional[QWidget]:
        """Get the main UI widget"""
        if self.ui_widget is None:
            try:
                from .visualizer_widget import VisualizerWidget
                self.ui_widget = VisualizerWidget(self)
            except Exception as e:
                logging.error(f"Failed to create widget: {e}")
                return None
        return self.ui_widget
    
    def cleanup(self):
        """Clean up plugin resources"""
        if self.ui_widget:
            self.ui_widget.deleteLater()
            self.ui_widget = None

# Plugin factory function
Plugin = ADOBMDVisualizer