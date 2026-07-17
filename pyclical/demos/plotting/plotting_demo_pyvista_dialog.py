#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# plotting_demo_pyvista_dialog.py: Interactive PySide6 dialog for PyVista plotting.
#
#    copyright            : (C) 2026 by Paul C. Leopardi
#
#    This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.

"""
An example of how to modify the PyVista visualization data via an interactive PySide6 dialog.
"""
import sys
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from pyvistaqt import QtInteractor

from plotting_demo_pyvista import (
    default_nbr_orbits,
    default_nbr_points,
    default_opacity,
    default_radius,
    default_resolution,
    default_scaling,
    default_segment_len,
    demo,
)


# pylint: disable=too-many-instance-attributes
class TorusDemoDialog(QMainWindow):
    """
    Interactive PySide6 dialog for PyClical PyVista torus demo.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Torus demo (PyVista)")
        self.resize(1280, 800)

        # Main widget and layout
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # Left panel for controls
        controls_panel = QWidget()
        controls_layout = QVBoxLayout(controls_panel)
        controls_panel.setFixedWidth(340)

        # Group 1: Plotting controls
        group_plotting = QGroupBox("Plotting controls")
        layout_plotting = QVBoxLayout(group_plotting)

        self.slider_orbits, self.label_orbits = self._create_slider(
            "Number of orbits:", 1, 16, default_nbr_orbits
        )
        layout_plotting.addLayout(self.label_orbits[0])
        layout_plotting.addWidget(self.slider_orbits)

        self.slider_res, self.label_res = self._create_slider(
            "Resolution:", 3, 32, default_resolution
        )
        layout_plotting.addLayout(self.label_res[0])
        layout_plotting.addWidget(self.slider_res)

        self.check_reciprocal = QCheckBox("Use reciprocal")
        self.check_reciprocal.setChecked(True)
        layout_plotting.addWidget(self.check_reciprocal)

        self.btn_plot = QPushButton("Draw a new plot")
        self.btn_plot.clicked.connect(self.new_plot)
        layout_plotting.addWidget(self.btn_plot)

        controls_layout.addWidget(group_plotting)

        # Group 2: Advanced controls
        group_adv = QGroupBox("Advanced controls")
        layout_adv = QVBoxLayout(group_adv)

        self.slider_points, self.label_points = self._create_slider(
            "Number of points:", 2, 32000, default_nbr_points
        )
        layout_adv.addLayout(self.label_points[0])
        layout_adv.addWidget(self.slider_points)

        self.slider_seglen, self.label_seglen = self._create_slider(
            "Segment length:", 1, 32000, default_segment_len
        )
        layout_adv.addLayout(self.label_seglen[0])
        layout_adv.addWidget(self.slider_seglen)

        self.slider_scaling, self.label_scaling = self._create_slider(
            "Scaling:", 1, 32000, default_scaling
        )
        layout_adv.addLayout(self.label_scaling[0])
        layout_adv.addWidget(self.slider_scaling)

        controls_layout.addWidget(group_adv)
        controls_layout.addStretch()

        main_layout.addWidget(controls_panel)

        # Right panel: embedded PyVista QtInteractor
        self.plotter = QtInteractor(self)
        main_layout.addWidget(self.plotter.interactor)

        # Initial plot
        self.new_plot()

    def _create_slider(self, label_text, min_val, max_val, default_val):
        h_layout = QHBoxLayout()
        lbl_title = QLabel(label_text)
        lbl_val = QLabel(str(default_val))
        h_layout.addWidget(lbl_title)
        h_layout.addStretch()
        h_layout.addWidget(lbl_val)

        slider = QSlider(Qt.Horizontal)
        slider.setMinimum(min_val)
        slider.setMaximum(max_val)
        slider.setValue(default_val)
        slider.valueChanged.connect(lambda v: lbl_val.setText(str(v)))

        return slider, (h_layout, lbl_val)

    def new_plot(self):
        demo(
            nbr_figures=1,
            nbr_orbits=self.slider_orbits.value(),
            nbr_points=self.slider_points.value(),
            scaling=self.slider_scaling.value(),
            segment_len=self.slider_seglen.value(),
            radius=default_radius,
            resolution=self.slider_res.value(),
            opacity=default_opacity,
            reciprocal=self.check_reciprocal.isChecked(),
            plotter=self.plotter,
        )
        self.plotter.reset_camera()

    def closeEvent(self, event):
        self.plotter.close()
        super().closeEvent(event)


def main():
    import os

    off_screen = (
        os.environ.get("PYVISTA_OFF_SCREEN", "").lower() in ("true", "1")
        or "--off-screen" in sys.argv
    )
    app = QApplication.instance() or QApplication(sys.argv)
    window = TorusDemoDialog()
    if off_screen:
        window.plotter.off_screen = True
        app.processEvents()
        window.close()
    else:
        window.show()
        sys.exit(app.exec())


if __name__ == "__main__":
    main()
