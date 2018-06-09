# -*- coding: utf-8 -*-
import dash
import dash_html_components as html
import dash_core_components as dcc
import tkinter as tk
from tkinter import filedialog

class App():
    def __init__(self):
        self.app = dash.Dash()
        self.folder = ""
        # self.root = tk.Tk()

        self.layout()

    def layout(self):
        self.app.layout = html.Div([
            html.Div(dcc.Input(id='input-box', type='text')),
            html.Button('Submit', id='button'),
            html.Div(id='output-container-button',
                     children='Enter a value and press submit')
        ])

        @self.app.callback(
            dash.dependencies.Output('output-container-button', 'children'),
            [dash.dependencies.Input('button', 'n_clicks')],
            [dash.dependencies.State('input-box', 'value')])
        def update_output(n_clicks, value):
            # return 'The input value was "{}" and the button has been clicked {} times'.format(
            #     value,
            #     n_clicks
            # )
            if n_clicks is not None:
                self.folder = filedialog.askdirectory()
                print(self.folder)


def main():

    tk.Tk().withdraw() # Close the root window
    in_path = filedialog.askdirectory()
    print(in_path)

dash_app = App()
if __name__ == '__main__':
    dash_app.app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
    dash_app.app.run_server(debug=True)