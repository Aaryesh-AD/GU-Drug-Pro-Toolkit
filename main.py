'''
GU Drug Pro Toolkit
Author: Aaryesh Deshpande

// Copyright (c) 2024 Aaryesh Deshpande
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so.

'''

# Module imports
import os, sys, webbrowser
import tkinter as tk
from tkinter import filedialog, Text, Scrollbar, Toplevel
from tkinter import LEFT, RIGHT, Y, BOTH, END, Menu
import customtkinter as ctk
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import pubchempy as pcp
import matplotlib.pyplot as plt
from functools import lru_cache
from tkinter import messagebox

# Adding specific paths to the system path for module access
from props_ptr.export_info import export_mol_pdf, export_mol_xls
from props_ptr.export_descp import compute_and_export_descriptors
from props_ptr.graphs_1 import plot_radar_chart, graph_prop_1
from props_ptr.name_chem import chem_name_info
from props_ptr.mol_prop import basic_mol_props
from props_ptr.mol_struct import calculate_molecular_parameters
from props_ptr.mol_exp import exp_info
from props_ptr.database_link import get_database_links
from props_ptr.sa_dl import calculate_sa_score
from props_ptr.struct_mdc import struct_exp
from props_ptr.drug_fil import rules_drg
from props_ptr.ESOL_prop import predict_sol
from props_ptr.pharmadmet import predict_for_all_models

'''
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
The code from here onwards defines a class named MoleculeAnalyzerApp
which inherits from ctk.CTk, a custom tkinter class for enhanced GUI elements.
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
'''

class MoleculeAnalyzerApp(ctk.CTk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # Calls the constructor of the parent class
        self.database_links = {}  # Initialize an empty dictionary for database links
        self.title("GU Drug-Pro Toolkit")  # Set the window title
        self.init_splash_screen()  # Initialize and display the splash screen
        self.splash.wait_window()  # Wait for the splash screen to close
        self.init_gui()  # Initialize the rest of the GUI components

        
    def init_splash_screen(self):
        # Method to create and display a splash screen upon application start
        self.splash = Toplevel()  # Creates a new top-level window
        self.splash.title("Loading...")  # Title for the splash screen
        splash_image_path = "asset/splash_400x.png"  # Path to the splash screen image
        pil_image = Image.open(splash_image_path)  # Open the image using PIL
        # Convert the PIL image to a customtkinter image
        splash_image = ctk.CTkImage(light_image=pil_image, size=(401,301))
        splash_label = ctk.CTkLabel(self.splash, text="", image=splash_image)
        splash_label.image = splash_image  # Keep a reference to avoid garbage collection
        splash_label.pack()  # Pack the label into the splash window
        # Set the geometry of the splash screen to center it on the screen
        width, height = 480, 360  # Dimensions of the splash screen
        screen_width = self.splash.winfo_screenwidth()  # Get the screen width
        screen_height = self.splash.winfo_screenheight()  # Get the screen height
        x = (screen_width / 2) - (width / 2)  # Calculate the x position
        y = (screen_height / 2) - (height / 2)  # Calculate the y position
        self.splash.geometry(f"{width}x{height}+{int(x)}+{int(y)}")  # Set the geometry
        self.splash.overrideredirect(True)  # Remove window decorations
        self.splash.after(2500, self.splash.destroy)  # Schedule the window to close after 2.5 seconds

    def about_init(self):
        if hasattr(self, 'about_window') and self.about_window.winfo_exists():
            self.about_window.lift()
            return
        self.about_window = ctk.CTkToplevel()  # Assign to self.about_window
        ico_path = "asset/icon.ico"  # Path to the window icon
        self.about_window.iconbitmap(ico_path)
        self.about_window.title("About")
        self.about_window.geometry("800x600")
        self.about_window.resizable(False, False)
        about_img = "asset/about_800x.png"
        pil_about = Image.open(about_img)
        about_image = ctk.CTkImage(light_image=pil_about, size=(801,601))
        about_label = ctk.CTkLabel(self.about_window, text="", image=about_image)
        about_label.image = about_image  # Keep a reference to avoid garbage collection
        about_label.pack()
    
    def export_desc(self):
        if hasattr(self, 'export_desc_win') and self.export_desc_win.winfo_exists():
            self.export_desc_win.lift()
            return    
        self.export_desc_win = ctk.CTkToplevel()
        ico_path = "asset/icon.ico"
        self.export_desc_win.iconbitmap(ico_path)
        self.export_desc_win.title("Export Descriptors")
        self.export_desc_win.geometry("350x300")
        self.export_desc_win.resizable(False, False)
        label = ctk.CTkLabel(self.export_desc_win, text="Export Descriptors Options", font=("Helvetica", 14))
        label.pack(pady=(5,0))
        # Checkboxes
        self.cb_morgan = ctk.CTkCheckBox(self.export_desc_win, text="2D Descriptors")
        self.cb_morgan.pack(pady=(10, 0))
        self.cb_ecfp = ctk.CTkCheckBox(self.export_desc_win, text="ECFP(X)")
        self.cb_ecfp.pack(pady=(5, 0))
        self.cb_maccs = ctk.CTkCheckBox(self.export_desc_win, text="MACCs")
        self.cb_maccs.pack(pady=(5, 0))
        self.cb_mqn = ctk.CTkCheckBox(self.export_desc_win, text="MQN")
        self.cb_mqn.pack(pady=(5, 0))
        # Dropdown for selecting CSV or Excel
        self.file_format_combobox = ctk.CTkComboBox(self.export_desc_win, values=["csv", "xlsx"])
        self.file_format_combobox.pack(pady=(20,5))

        def export_action():
            query = self.smiles_entry.get()
            smiles = self.retrieve_smiles(query)
            if not smiles:
                messagebox.showerror("Error", "Please enter a valid SMILES string.")
                return
            
            selected_format = self.file_format_combobox.get()
            
            # Open a dialog box for the user to select a directory
            export_directory = filedialog.askdirectory()
            if not export_directory:
                messagebox.showerror("Error", "No directory selected.")
                return
            
            # Check for descriptor selection
            options = {
                '2D': self.cb_morgan.get(),
                'ECFP': self.cb_ecfp.get(),
                'MACCS': self.cb_maccs.get(),
                'MQN': self.cb_mqn.get()
            }
            
            if not any(options.values()):  # If no options are selected
                messagebox.showerror("Error", "Please select at least one descriptor type.")
                return
            
            # Pass the directory to the export function
            try:
                compute_and_export_descriptors(smiles, selected_format, options, export_directory)
                messagebox.showinfo("Success", f"Data successfully exported as {selected_format}.")
            except Exception as e:
                messagebox.showerror("Export Failed", str(e))
            finally:
                self.export_desc_win.destroy() 
        export_button = ctk.CTkButton(self.export_desc_win, text="Export", command=export_action)
        export_button.pack(pady=(10, 0), padx=10)
        cancel_button = ctk.CTkButton(self.export_desc_win, text="Cancel",fg_color="#9b3521", text_color="white", hover_color="#ce472c", command=self.export_desc_win.destroy)
        cancel_button.pack(pady=(5, 0), padx = 10) 
    
    def open_export_window(self):
        if hasattr(self, 'export_win') and self.export_win.winfo_exists():
            self.export_win.lift()
            return 
        self.export_win = ctk.CTkToplevel()
        ico_path = "asset/icon.ico"
        self.export_win.iconbitmap(ico_path)
        self.export_win.title("Export Result Data")
        self.export_win.geometry("350x150")
        self.export_win.resizable(False, False)

        label = ctk.CTkLabel(self.export_win, text="Select Export Format", font=("Helvetica", 14))
        label.pack(pady=(10, 5))

        # Dropdown for selecting the export format
        self.export_format_combobox = ctk.CTkComboBox(self.export_win, values=["Excel", "PDF"])
        self.export_format_combobox.pack(pady=(0, 20))

        # Export button
        export_button = ctk.CTkButton(self.export_win, text="Export", command=self.export_data)
        export_button.pack(side=tk.LEFT, padx=20, pady=10)

        # Cancel button
        cancel_button = ctk.CTkButton(self.export_win, text="Cancel", fg_color="#9b3521", text_color="white", hover_color="#ce472c", command=self.export_win.destroy)
        cancel_button.pack(side=tk.RIGHT, padx=20, pady=10)
        
    def export_data(self):
        query = self.smiles_entry.get()
        smiles = self.retrieve_smiles(query) # Retrieve the SMILES string from the entry
        if not smiles:
            messagebox.showerror("Error", "Please enter a valid SMILES string.")
            return
        
        export_format = self.export_format_combobox.get()
        if export_format not in ["Excel", "PDF"]:
            messagebox.showerror("Error", "Please select a valid export format.")
            return

        # Set up the file types for the save dialog based on the selected format
        file_types = ("Excel files", "*.xlsx") if export_format == "Excel" else ("PDF files", "*.pdf")
        filename = filedialog.asksaveasfilename(defaultextension=file_types[1], filetypes=[file_types], initialfile = f"{smiles}_result data" )
        if not filename:
            messagebox.showerror("Error", "No file name specified.")
            return
        
        try:
            if export_format == "Excel":
                export_mol_xls(smiles, filename)  # Call function to export to Excel
            elif export_format == "PDF":
                export_mol_pdf(smiles, filename)  # Call function to export to PDF
            messagebox.showinfo("Success", f"Data successfully exported as {export_format}.")
        except Exception as e:
            messagebox.showerror("Export Failed", str(e))
        finally:
            self.export_win.destroy()  # Close the export window after operation  
        
    def on_close(self):
        plt.close('all')  # Close all matplotlib figures
        # If you have a specific FigureCanvasTkAgg widget, destroy it here
        # e.g., self.canvas.get_tk_widget().destroy()
        self.destroy() 
        
    def init_gui(self):
        # Method to initialize the main GUI components after the splash screen
        self.update_idletasks()  # Update the GUI state
        self.configure_gui()  # Configure the main window's appearance and geometry
        self.create_menu_bar()  # Create the menu bar with its items
        self.create_upload_frame()  # Create the frame for uploading molecules
        self.create_display_frame()  # Create the frame for displaying molecule information
        self.create_tab_view()  # Create the tabbed view for different information categories
        self.pack_gui()  # Pack the main GUI components into the window
        self.protocol("WM_DELETE_WINDOW", self.on_close)
        
    def configure_gui(self):
        # Set the appearance mode and window geometry
        ctk.set_appearance_mode("Light")  # Set the theme
        width = self.winfo_screenwidth()  # Get the width of the screen
        height = self.winfo_screenheight()  # Get the height of the screen
        self.geometry(f"{width}x{height}")  # Set the window size to full screen
        ico_path = "asset/icon.ico"  # Path to the window icon
        self.iconbitmap(ico_path)  # Set the window icon
    
    def open_doc(self):
        import webbrowser
        pdf_path = 'asset/manual.pdf'
        absolute_path = os.path.abspath(pdf_path)
        url = f'file://{absolute_path}'
        webbrowser.open(url)

    def create_menu_bar(self):
        # Method to create a menu bar at the top of the main window
        menu_bar = Menu(self)  # Create a new menu bar
        self.config(menu=menu_bar)  # Attach the menu bar to the main window
        # Create and add the drop-down menus to the menu bar
        for label in ["File", "Export", "Help"]:
            menu = Menu(menu_bar, tearoff=0, font=("",11))  # Create a new drop-down menu
            if label == 'File':
                menu.add_command(label= "Clear and Refresh", command= self.clear_input)# Add a dummy command to the menu
                menu.add_separator()
                menu.add_command(label= "Exit", command = self.on_close)
            if label == "Export":
                menu.add_command(label= "Export Descriptors", command= self.export_desc)
                menu.add_separator()
                menu.add_command(label= "Export Result Data", command=self.open_export_window)
            if label == "Help":
                menu.add_command(label="About", command=self.about_init)
                menu.add_separator()
                menu.add_command(label="Documentation", command = self.open_doc)
            menu_bar.add_cascade(label=label, menu=menu)  # Add the menu to the menu bar
            
    def create_upload_frame(self):
        # Create a frame for uploading molecules or entering SMILES strings
        self.row_frame = ctk.CTkFrame(self)  # Create a new frame
        self.row_frame.grid_columnconfigure(2, weight=1)  # Configure the grid layout
        # Create and configure the upload button
        self.upload_button = ctk.CTkButton(
            self.row_frame, text="Upload Molecule", command=self.upload_molecule,
            border_width=2, border_color="#777", fg_color="#dbdbdb",
            text_color="#777", hover_color="#b5b5b5"
        )
        self.upload_button.grid(row=0, column=0, padx=10)  # Place the button in the grid
        # Create and configure a label indicating an alternative input method
        option_title = ctk.CTkLabel(self.row_frame, text="or")
        option_title.grid(row=0, column=1, padx=10)  # Place the label in the grid
        # Create and configure the SMILES entry widget
        self.smiles_entry = ctk.CTkEntry(
            self.row_frame,
            placeholder_text="Enter SMILES / IUPAC Name / Chemical Name / Molecular Formula / InChi......... "
        )
        self.smiles_entry.grid(row=0, column=2, sticky='ew', padx=10)  # Place the entry widget in the grid
        # Create and configure the calculate button
        self.calc_button = ctk.CTkButton(self.row_frame, text="Calculate", command=self.show_molecule)
        self.calc_button.grid(row=0, column=3, padx=10)  # Place the button in the grid
        # Create and configure the clear button
        self.clear_button = ctk.CTkButton(
            self.row_frame, text="Clear", command=self.clear_input,
            fg_color="#9b3521", text_color="white", hover_color="#ce472c"
        )
        self.clear_button.grid(row=0, column=4, padx=10)  # Place the button in the grid
        
    def create_display_frame(self):
        # Create a frame for displaying the molecule image and properties
        self.display_frame = ctk.CTkFrame(self)  # Create a new frame
        self.display_frame.grid_columnconfigure(1, weight=3)  # Configure the grid layout
        # Create and configure the label for displaying the molecule image
        self.display_label = ctk.CTkLabel(self.display_frame, text="", width=180, height=180)
        self.display_label.grid(row=0, column=0, padx=10, pady=10)  # Place the label in the grid
        # Create a frame for displaying molecular properties
        self.prop_frame = ctk.CTkFrame(self.display_frame)
        self.prop_frame.grid(row=0, column=1, sticky="nsew")  # Place the frame in the grid
        # Create and configure the title label for molecular properties
        self.prop_name_title = ctk.CTkLabel(
            self.prop_frame, text='', font=("Helvetica", 16, "bold"), justify=LEFT, anchor="w"
        )
        self.prop_name_title.pack(pady=(20,0), padx=30, fill='x')  # Pack the label into the frame
        # Create and configure the label for displaying the properties themselves
        self.name_prop_label = ctk.CTkLabel(
            self.prop_frame, text='', font=("Helvetica", 12), justify=LEFT, anchor="w", wraplength=1100
        )
        self.name_prop_label.pack(pady=(5,20), padx=30, fill='x')  # Pack the label into the frame
        # Container for database link buttons
        self.db_link_container = ctk.CTkFrame(self.prop_frame)
        self.db_link_container.pack(pady=(5,20), padx=30, fill='x', expand=True)
        
    def create_tab_view(self):
        # Create a tabbed view for displaying various categories of information
        self.tab_view_frame = ctk.CTkFrame(self)  # Create a new frame for the tabbed view
        self.tab_view = ctk.CTkTabview(self.tab_view_frame)  # Create the tabbed view
        # Add tabs for different categories of information
        for tab_name in ["Structural Info", "Record Data", "Physicochemical Info", "Pharmacokinetics", "Druglikeness", "Stats"]:
            tab = self.tab_view.add(tab_name)  # Add a new tab
            if tab_name == "Stats":
                self.create_stat_view(tab)
            else:
                # Create a text widget for displaying information in the tab
                text_widget = Text(tab, wrap='word', font=("Helvetica", 11), padx=10, pady=10)
                text_widget.pack(side=LEFT, fill=BOTH, expand=True)  # Pack the text widget into the tab
                # Create a scrollbar for the text widget
                scrollbar = Scrollbar(tab, command=text_widget.yview)
                scrollbar.pack(side=RIGHT, fill=Y)  # Pack the scrollbar into the tab
                text_widget.config(yscrollcommand=scrollbar.set)  # Link the scrollbar to the text widget
                def ignore_keypress(event):
                    return "break"  # Prevent text input in the widget
                text_widget.bind("<KeyPress>", ignore_keypress)  # Bind the key press event
                # Dynamically set an attribute for easy access to the text widget later
                setattr(self, f"text_{tab_name.replace(' ', '_').lower()}", text_widget)
            
    def clear_tab_content(self, tab):
        if tab.winfo_exists():
            for widget in tab.winfo_children():
                widget.destroy()  
                
    def create_stat_view(self, parent_tab, smiles=None):
        # Initialize the statistics view and tabs only once
        if not hasattr(self, 'stat_view_initialized') or not self.stat_view_initialized:
            self.stat_view_frame = ctk.CTkFrame(parent_tab)
            self.stat_view = ctk.CTkTabview(self.stat_view_frame)
            # Add tabs for the spider web chart and other graphs
            self.spiderweb_tab = self.stat_view.add("Bioavailability Web")
            self.stat_view_initialized = True
            self.stat_view_frame.pack(fill='both', expand=True)
            self.stat_view.pack(fill='both', expand=True, padx=10, pady=10)

        if smiles:
            # If smiles is provided, plot the graph in the appropriate tab
            compound_values = graph_prop_1(smiles)
            # Update the Spiderweb chart tab with the new content
            self.clear_tab_content(self.spiderweb_tab)  # Clear previous content if any
            canvas, fig = plot_radar_chart(self.spiderweb_tab, compound_values)
            canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
            # No     further action needed for 'self.graph_tab' unless you plan to plot something different there

    def update_database_links(self):
        # Method to update the database links displayed in the GUI
        for widget in self.db_link_container.winfo_children():
            widget.destroy()  # Remove existing widgets
        # Add new buttons for each database link
        button_frame = ctk.CTkFrame(self.db_link_container)
        button_frame.pack(pady=10, expand=True)  # Pack the frame for buttons
        for db_name, db_url in self.database_links.items():
            button = ctk.CTkButton(
                button_frame, text=db_name, command=lambda url=db_url: webbrowser.open(url),
                border_width=2, border_color="#444444", fg_color="#dbdbdb",
                text_color="#444444", hover_color="#b5b5b5"
            )
            button.pack(side=tk.LEFT, padx=5, anchor='center')  # Pack each button

    def pack_gui(self):
        # Method to pack the main GUI components into the window
        self.row_frame.pack(pady=(15, 10), padx=20, fill='x')  # Pack the upload frame
        self.display_frame.pack(pady=(10, 10), padx=20, fill='both', expand=False)  # Pack the display frame
        self.tab_view_frame.pack(pady=(0, 10), padx=20, fill='both', expand=True)  # Pack the tab view frame
        self.tab_view.pack(expand=True, fill="both", padx=10, pady=10)  # Pack the tab view
        
    def clear_input(self):
        # Method to reset the input fields and clear displayed information
        self.smiles_entry.delete(0, 'end')  # Clear the SMILES entry widget
        self.smiles_entry.configure(state="normal")  # Enable the SMILES entry widget
        self.upload_button.configure(state="normal")  # Enable the upload button
        self.display_label.configure(image="")  # Clear the molecule image
        self.display_label.image = None  # Clear the reference to the image to avoid memory leaks
        self.display_label.configure(text="")  # Clear any text in the display label
        self.name_prop_label.configure(text="")  # Clear the properties text
        self.prop_name_title.configure(text='')  # Clear the properties title
        # Clear text widgets in each information category tab
        self.text_structural_info.delete(1.0, tk.END)
        self.text_physicochemical_info.delete(1.0, tk.END)
        self.text_pharmacokinetics.delete(1.0, tk.END)
        self.text_druglikeness.delete(1.0, tk.END)
        self.text_record_data.delete(1.0, tk.END)
        self.database_links = {}  # Reset the database links dictionary
        self.update_database_links()  # Update the UI to reflect the cleared database links
        self.smiles_entry.delete(0, 'end')  # Clear any uncleared SMILES from entry widget
        self.smiles_entry.configure(state="normal")  # Enable the SMILES entry widget
        self.clear_tab_content(self.spiderweb_tab)
        
    @lru_cache(maxsize=128)
    def retrieve_smiles(self, query):
        """Retrieve SMILES string for a given query using PubChem."""
        # Try to interpret the query as a SMILES string directly
        if Chem.MolFromSmiles(query) is not None:
            return query  # Return the query if it's a valid SMILES string
        # Try to fetch the SMILES string from PubChem if the query is not a valid SMILES string
        try:
            compound = pcp.get_compounds(query, 'name')
            if compound:
                return compound[0].isomeric_smiles  # Return the isomeric SMILES string
        except pcp.PubChemHTTPError:
            tk.messagebox.showerror("Network Error", "Check your network connection and try again.")
            return None
        except Exception as e:
            tk.messagebox.showerror("Error", str(e))
            return None
        tk.messagebox.showerror("Error", "Compound not found or invalid query.")
        return None  # Return None if no valid SMILES string could be retrieved
        
    def upload_molecule(self):
        """Handle molecule upload from a file."""
        self.clear_input()  # Reset input fields and clear displayed information
        file_path = filedialog.askopenfilename(filetypes=[("MOL files", "*.mol"), ("SDF files", "*.sdf")])
        if file_path:
            # Load the molecule from the specified file
            if file_path.endswith('.mol'):
                mol = Chem.MolFromMolFile(file_path)
            elif file_path.endswith('.sdf'):
                suppl = Chem.SDMolSupplier(file_path)
                mol = next(suppl)
            else:
                return
            if mol is not None:
                # Convert the molecule to a SMILES string
                smiles = Chem.MolToSmiles(mol)
                self.smiles_entry.delete(0, 'end')  # Clear the SMILES entry widget
                self.smiles_entry.insert(0, smiles)  # Insert the SMILES string into the entry widget
                self.smiles_entry.configure(state="disabled")
                query = self.smiles_entry.get()
                if query:
                    smiles = self.retrieve_smiles(query)# Disable the SMILES entry widget
                if smiles:  
                    self.show_molecule(smiles)  # Display the molecule based on the SMILES string
            else:
                self.clear_input()  # Reset input fields and clear displayed information if loading fails
                
    def show_molecule(self, smiles=None):
        """Display information about the molecule specified by a SMILES string."""
        if smiles is None:
            # Retrieve the SMILES string from the entry widget if not provided
            query = self.smiles_entry.get()
            if query:
                smiles = self.retrieve_smiles(query)
        if smiles:
            molecule = Chem.MolFromSmiles(smiles)
            if self.smiles_entry.cget('state') == "normal":
                self.upload_button.configure(state="disabled")  # Disable the upload button if input is valid
            if molecule:
                # Render the molecule to an image and display it
                self.display_label.configure(text="")
                img = Draw.MolToImage(molecule, size=(300, 200))
                img.save("temp_img/molecule.png")  # Save the rendered image to a file
                img_mol = Image.open("temp_img/molecule.png")  # Load the image file
                photo = ctk.CTkImage(light_image=img_mol, size=(300,200))  # Convert to a customtkinter image
                self.display_label.configure(image=photo)  # Display the molecule image
                self.display_label.image = photo  # Keep a reference to avoid garbage collection
                # Retrieve and display various properties and information about the molecule
                params_name = chem_name_info(smiles)
                mol_prop_func = basic_mol_props(smiles)
                mol_struct_func = calculate_molecular_parameters(smiles)
                properties_text = "\n\n".join(f"{key}: {value}" for key, value in params_name.items())
                mol_prop_name = "\n\n".join(f"{key}: {value}" for key, value in mol_prop_func.items())
                mol_struct_name = "\n\n".join(f"{key}:\t{value}" for key, value in mol_struct_func.items())
                struct_name = mol_prop_name + "\n\n" + mol_struct_name
                sa_func = calculate_sa_score(smiles)
                sa_name = "\n\n".join(f"{key}: {value}" for key, value in sa_func.items())
                drg_struc_func = struct_exp(smiles)
                drg_struc_name = "\n\n".join(f"{key}: {value}" for key, value in drg_struc_func.items())
                drg_lik_func = rules_drg(smiles)
                drg_lik_name = "\n\n".join(f"{key}: {value}" for key, value in drg_lik_func.items())
                drg_prop = sa_name + "\n\n" + drg_struc_name + "\n\n" + drg_lik_name
                phys_func = predict_sol(smiles)
                phys_name = "\n\n".join(f"{key}: {value}" for key, value in phys_func.items())
                exp_name = exp_info(smiles)
                exp_func = "\n\n".join(f"{key}: {value}" for key, value in exp_name.items())
                adme_name = predict_for_all_models(smiles)
                adme_func = "\n\n".join(f"{key}: {value}" for key, value in adme_name.items())
                # Update the database links based on the molecule
                self.database_links = get_database_links(smiles)
                self.update_database_links()
                # Display the names and identifiers
                self.name_prop_label.configure(text=properties_text)
                self.prop_name_title.configure(text='Names and Identifiers')
                # Update the text widgets in each tab with the relevant information
                self.text_structural_info.delete(1.0, END)
                self.text_structural_info.insert(END, struct_name)
                self.text_physicochemical_info.delete(1.0, END)
                self.text_physicochemical_info.insert(END, phys_name)
                self.text_pharmacokinetics.delete(1.0, END)
                self.text_pharmacokinetics.insert(END, adme_func)  # Placeholder, replace as needed
                self.text_druglikeness.delete(1.0, END)
                self.text_druglikeness.insert(END, drg_prop)
                self.text_record_data.delete(1.0, END)
                self.text_record_data.insert(END, exp_func)
                self.create_stat_view(self.stat_view_frame, smiles)



            else:
                # Handle invalid input by clearing the display and showing an error message
                self.display_label.configure(image="")
                self.display_label.configure(text="Invalid SMILES / Mol file")
                tk.messagebox.showerror("Error", "Invalid SMILES / Mol file")
                self.name_prop_label.configure(text="")
                self.prop_name_title.configure(text='')
        else:
            # Handle missing or invalid input by clearing the display and prompting for valid input
            self.display_label.configure(image="")
            self.display_label.configure(text="Enter SMILES, IUPAC Name, Chemical Name, Molecular Formula, InChi.... / Upload Mol file")
            self.name_prop_label.configure(text="")
            self.prop_name_title.configure(text='')

if __name__ == "__main__":
    app = MoleculeAnalyzerApp()  # Create an instance of the application
    app.mainloop()  # Start the application's main event loop
