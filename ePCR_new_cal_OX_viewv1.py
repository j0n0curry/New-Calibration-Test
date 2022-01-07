

import streamlit as st
import os
from io import BytesIO
from io import StringIO
from io import TextIOWrapper
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math
 
 
st.set_page_config(layout="wide")

 


version = 'New Cal - Oxford - v1.1 - alpha'

def main():
    
    
    #Set up main page of application / Header info / data collection / file selection / remove files / Reset
    
    #Main landing page greating / info
    

    st.title('ePCR analysis tool ' +str(version))

    st.subheader("Upload Araya csv file directly for processing - either drag and drop or click to select files on your local machine")
    
    
    
  
    
    #Select main function iinstatiat ArayaManager before for loop otherwise reads and closes files

    data_manager = ArayaManager()

    #upload files - user select files to upload to server from local system / or other system
    # files are bytes like therefore you read and open and close the file during running.
    #unexpected behaviour like loss of data / loss of functions - especially int and string replacements. 
    #uploaded files is a list of file like objects - not diretly readable by Pandas. 
    with st.form("upload", clear_on_submit=True):
        uploaded_files = st.file_uploader("Choose a CSV file", type=['csv'], accept_multiple_files=True, key = 'key', help = 'select files by highlighting them then press upload file to start process' )
        submitted = st.form_submit_button("UPLOAD!",  help = 'upload files to start analysis!')

 
   
        
    #loop through uploaded files to break open the list and pass through data IO stream 
    for uploaded_file in uploaded_files:
        if submitted and uploaded_file is not None: # stops the error of 'no files and start function if there are uploaed files
            file_names = (uploaded_file.name)# to add to Araya Manager.
            data_manager.get_run_name(file_names) # get the Run_ID. access ArayaManager function by passing file name
            data_manager.get_date_time(file_names)
            data_manager.concatenate_dataframes(uploaded_file) # concat dataframes.
            
            
    comp = data_manager.group_df #assign variable - contactenated dataframes - csv loading direct from path doesn't require this. 
    print(comp)
    
    # start onwards with processing only if dataframe us greater than 0 
    if len(comp) > 0:
        #conversion of string to float in ArayaManager needed - test adding in a coersion function - as per below. 
        #coerce mixed float / int nummbers from somewhere. Add to Araya Manager a method coerce_numeric.
        comp['FAM_RFU'] = comp['FAM_RFU'].astype('float').astype('Int32')
        comp['VIC_RFU'] = comp['VIC_RFU'].astype('float').astype('Int32')
        comp['ROX_RFU'] = comp['ROX_RFU'].astype('float').astype('Int32')
        #will remove this function - to the bottom - allow files to be processed and user defined and downloaed as whole set and analysed sets. 
        
        #process file attributes in to parameters for QC. Essential information. 
        comp['Well'] = comp['Row_ID']+comp['Col_ID']
        comp['norm_RNaseP'] =  comp['VIC_RFU'].abs() / comp['ROX_RFU']
        comp['norm_N_Cov'] =  comp["FAM_RFU"]  / comp['ROX_RFU']
        comp.index.names=['order']
        comp.reset_index(inplace = True)
        comp['date_time'] = pd.to_datetime(comp['date_time'], format='%Y%m%d%H%M%S')
        comp[['date', 'time']] = comp['date_time'].astype(str).str.split(' ', 1, expand=True)
        print(comp.norm_RNaseP)
        comp.to_csv('test_out.csv')
        
        
        
        controls = {'P19': 'A1500', 'O19' : 'A1500', 'O20': 'A1500',
                        'P21': 'NEG', 'O21': 'NEG', 'O22': 'NEG',
                        'O23':'S06', 'P23':'S06', 'O24': 'S06'}
           
        comp['control'] = comp['Well'].map(controls).fillna('paitent')
      
    def valid_array(df):
        vdf = df.groupby('Run_ID')['ROX_RFU'].mean().reset_index()
        vdf = vdf[vdf['ROX_RFU'] <= 1000] 
        vdf =  vdf['Run_ID'].to_list()
        df = df[~df.Run_ID.isin(vdf)]
        st.write('Arrays Exluded', str(vdf))
        return(df)    
        
    comp=valid_array(comp)    
        
    def scoring(row):
    
        if row['norm_N_Cov'] < 0.5 and row['norm_RNaseP'] > 0.2:
            return('Negative Patient')
        elif row['norm_N_Cov'] > 0.5 and row['norm_N_Cov'] <= 1.0 and row['norm_RNaseP'] >0.2:
            return('PLOD')
        elif row['norm_N_Cov'] > 1.0 and row['norm_RNaseP'] >= 0.2:
            return('N_Cov Paitent Positive')
        elif row['norm_N_Cov'] > 1.0 and row['norm_RNaseP']<= 0.2:
            return('Control_N_Cov')
        elif row['norm_N_Cov'] <= 0.5 and row['norm_RNaseP'] <=0.2:
            return('No_Call')
        elif row['norm_N_Cov'] > 0.5 and row['norm_N_Cov'] <= 1.0 and row['norm_RNaseP'] <0.2:
            return'CTRL_PLOD'
        else:
            return('missing')
        
    def void_detect_neg(row):
        if row['norm_N_Cov'] > 1.0:
            return('Detected')
        elif row['norm_N_Cov'] > 0.5 and row['norm_N_Cov'] <= 1.0:
            return('PLOD')
        else:
            return('negative')
    
         
    
    
    comp['Result'] = comp.apply(lambda row: scoring(row), axis = 1)   
    comp['Detection'] = comp.apply(lambda row: void_detect_neg(row), axis = 1)   
    
    quad_384 = {'A1':'QUAD 1','A2':'QUAD 2','A3':'QUAD 1','A4':'QUAD 2','A5':'QUAD 1','A6':'QUAD 2','A7':'QUAD 1','A8':'QUAD 2','A9':'QUAD 1','A10':'QUAD 2','A11':'QUAD 1','A12':'QUAD 2','A13':'QUAD 1','A14':'QUAD 2','A15':'QUAD 1','A16':'QUAD 2','A17':'QUAD 1','A18':'QUAD 2','A19':'QUAD 1','A20':'QUAD 2','A21':'QUAD 1','A22':'QUAD 2','A23':'QUAD 1','A24':'QUAD 2',
                'B1':'QUAD 3','B2':'QUAD 4','B3':'QUAD 3','B4':'QUAD 4','B5':'QUAD 3','B6':'QUAD 4','B7':'QUAD 3','B8':'QUAD 4','B9':'QUAD 3','B10':'QUAD 4','B11':'QUAD 3','B12':'QUAD 4','B13':'QUAD 3','B14':'QUAD 4','B15':'QUAD 3','B16':'QUAD 4','B17':'QUAD 3','B18':'QUAD 4','B19':'QUAD 3','B20':'QUAD 4','B21':'QUAD 3','B22':'QUAD 4','B23':'QUAD 3','B24':'QUAD 4',
                'C1':'QUAD 1','C2':'QUAD 2','C3':'QUAD 1','C4':'QUAD 2','C5':'QUAD 1','C6':'QUAD 2','C7':'QUAD 1','C8':'QUAD 2','C9':'QUAD 1','C10':'QUAD 2','C11':'QUAD 1','C12':'QUAD 2','C13':'QUAD 1','C14':'QUAD 2','C15':'QUAD 1','C16':'QUAD 2','C17':'QUAD 1','C18':'QUAD 2','C19':'QUAD 1','C20':'QUAD 2','C21':'QUAD 1','C22':'QUAD 2','C23':'QUAD 1','C24':'QUAD 2',
                'D1':'QUAD 3','D2':'QUAD 4','D3':'QUAD 3','D4':'QUAD 4','D5':'QUAD 3','D6':'QUAD 4','D7':'QUAD 3','D8':'QUAD 4','D9':'QUAD 3','D10':'QUAD 4','D11':'QUAD 3','D12':'QUAD 4','D13':'QUAD 3','D14':'QUAD 4','D15':'QUAD 3','D16':'QUAD 4','D17':'QUAD 3','D18':'QUAD 4','D19':'QUAD 3','D20':'QUAD 4','D21':'QUAD 3','D22':'QUAD 4','D23':'QUAD 3','D24':'QUAD 4',
                'E1':'QUAD 1','E2':'QUAD 2','E3':'QUAD 1','E4':'QUAD 2','E5':'QUAD 1','E6':'QUAD 2','E7':'QUAD 1','E8':'QUAD 2','E9':'QUAD 1','E10':'QUAD 2','E11':'QUAD 1','E12':'QUAD 2','E13':'QUAD 1','E14':'QUAD 2','E15':'QUAD 1','E16':'QUAD 2','E17':'QUAD 1','E18':'QUAD 2','E19':'QUAD 1','E20':'QUAD 2','E21':'QUAD 1','E22':'QUAD 2','E23':'QUAD 1','E24':'QUAD 2',
                'F1':'QUAD 3','F2':'QUAD 4','F3':'QUAD 3','F4':'QUAD 4','F5':'QUAD 3','F6':'QUAD 4','F7':'QUAD 3','F8':'QUAD 4','F9':'QUAD 3','F10':'QUAD 4','F11':'QUAD 3','F12':'QUAD 4','F13':'QUAD 3','F14':'QUAD 4','F15':'QUAD 3','F16':'QUAD 4','F17':'QUAD 3','F18':'QUAD 4','F19':'QUAD 3','F20':'QUAD 4','F21':'QUAD 3','F22':'QUAD 4','F23':'QUAD 3','F24':'QUAD 4',
                'G1':'QUAD 1','G2':'QUAD 2','G3':'QUAD 1','G4':'QUAD 2','G5':'QUAD 1','G6':'QUAD 2','G7':'QUAD 1','G8':'QUAD 2','G9':'QUAD 1','G10':'QUAD 2','G11':'QUAD 1','G12':'QUAD 2','G13':'QUAD 1','G14':'QUAD 2','G15':'QUAD 1','G16':'QUAD 2','G17':'QUAD 1','G18':'QUAD 2','G19':'QUAD 1','G20':'QUAD 2','G21':'QUAD 1','G22':'QUAD 2','G23':'QUAD 1','G24':'QUAD 2',
                'H1':'QUAD 3','H2':'QUAD 4','H3':'QUAD 3','H4':'QUAD 4','H5':'QUAD 3','H6':'QUAD 4','H7':'QUAD 3','H8':'QUAD 4','H9':'QUAD 3','H10':'QUAD 4','H11':'QUAD 3','H12':'QUAD 4','H13':'QUAD 3','H14':'QUAD 4','H15':'QUAD 3','H16':'QUAD 4','H17':'QUAD 3','H18':'QUAD 4','H19':'QUAD 3','H20':'QUAD 4','H21':'QUAD 3','H22':'QUAD 4','H23':'QUAD 3','H24':'QUAD 4',
                'I1':'QUAD 1','I2':'QUAD 2','I3':'QUAD 1','I4':'QUAD 2','I5':'QUAD 1','I6':'QUAD 2','I7':'QUAD 1','I8':'QUAD 2','I9':'QUAD 1','I10':'QUAD 2','I11':'QUAD 1','I12':'QUAD 2','I13':'QUAD 1','I14':'QUAD 2','I15':'QUAD 1','I16':'QUAD 2','I17':'QUAD 1','I18':'QUAD 2','I19':'QUAD 1','I20':'QUAD 2','I21':'QUAD 1','I22':'QUAD 2','I23':'QUAD 1','I24':'QUAD 2',
                'J1':'QUAD 3','J2':'QUAD 4','J3':'QUAD 3','J4':'QUAD 4','J5':'QUAD 3','J6':'QUAD 4','J7':'QUAD 3','J8':'QUAD 4','J9':'QUAD 3','J10':'QUAD 4','J11':'QUAD 3','J12':'QUAD 4','J13':'QUAD 3','J14':'QUAD 4','J15':'QUAD 3','J16':'QUAD 4','J17':'QUAD 3','J18':'QUAD 4','J19':'QUAD 3','J20':'QUAD 4','J21':'QUAD 3','J22':'QUAD 4','J23':'QUAD 3','J24':'QUAD 4',
                'K1':'QUAD 1','K2':'QUAD 2','K3':'QUAD 1','K4':'QUAD 2','K5':'QUAD 1','K6':'QUAD 2','K7':'QUAD 1','K8':'QUAD 2','K9':'QUAD 1','K10':'QUAD 2','K11':'QUAD 1','K12':'QUAD 2','K13':'QUAD 1','K14':'QUAD 2','K15':'QUAD 1','K16':'QUAD 2','K17':'QUAD 1','K18':'QUAD 2','K19':'QUAD 1','K20':'QUAD 2','K21':'QUAD 1','K22':'QUAD 2','K23':'QUAD 1','K24':'QUAD 2',
                'L1':'QUAD 3','L2':'QUAD 4','L3':'QUAD 3','L4':'QUAD 4','L5':'QUAD 3','L6':'QUAD 4','L7':'QUAD 3','L8':'QUAD 4','L9':'QUAD 3','L10':'QUAD 4','L11':'QUAD 3','L12':'QUAD 4','L13':'QUAD 3','L14':'QUAD 4','L15':'QUAD 3','L16':'QUAD 4','L17':'QUAD 3','L18':'QUAD 4','L19':'QUAD 3','L20':'QUAD 4','L21':'QUAD 3','L22':'QUAD 4','L23':'QUAD 3','L24':'QUAD 4',
                'M1':'QUAD 1','M2':'QUAD 2','M3':'QUAD 1','M4':'QUAD 2','M5':'QUAD 1','M6':'QUAD 2','M7':'QUAD 1','M8':'QUAD 2','M9':'QUAD 1','M10':'QUAD 2','M11':'QUAD 1','M12':'QUAD 2','M13':'QUAD 1','M14':'QUAD 2','M15':'QUAD 1','M16':'QUAD 2','M17':'QUAD 1','M18':'QUAD 2','M19':'QUAD 1','M20':'QUAD 2','M21':'QUAD 1','M22':'QUAD 2','M23':'QUAD 1','M24':'QUAD 2',
                'N1':'QUAD 3','N2':'QUAD 4','N3':'QUAD 3','N4':'QUAD 4','N5':'QUAD 3','N6':'QUAD 4','N7':'QUAD 3','N8':'QUAD 4','N9':'QUAD 3','N10':'QUAD 4','N11':'QUAD 3','N12':'QUAD 4','N13':'QUAD 3','N14':'QUAD 4','N15':'QUAD 3','N16':'QUAD 4','N17':'QUAD 3','N18':'QUAD 4','N19':'QUAD 3','N20':'QUAD 4','N21':'QUAD 3','N22':'QUAD 4','N23':'QUAD 3','N24':'QUAD 4',
                'O1':'QUAD 1','O2':'QUAD 2','O3':'QUAD 1','O4':'QUAD 2','O5':'QUAD 1','O6':'QUAD 2','O7':'QUAD 1','O8':'QUAD 2','O9':'QUAD 1','O10':'QUAD 2','O11':'QUAD 1','O12':'QUAD 2','O13':'QUAD 1','O14':'QUAD 2','O15':'QUAD 1','O16':'QUAD 2','O17':'QUAD 1','O18':'QUAD 2','O19':'QUAD 1','O20':'QUAD 2','O21':'QUAD 1','O22':'QUAD 2','O23':'QUAD 1','O24':'QUAD 2',
                'P1':'QUAD 3','P2':'QUAD 4','P3':'QUAD 3','P4':'QUAD 4','P5':'QUAD 3','P6':'QUAD 4','P7':'QUAD 3','P8':'QUAD 4','P9':'QUAD 3','P10':'QUAD 4','P11':'QUAD 3','P12':'QUAD 4','P13':'QUAD 3','P14':'QUAD 4','P15':'QUAD 3','P16':'QUAD 4','P17':'QUAD 3','P18':'QUAD 4','P19':'QUAD 3','P20':'QUAD 4','P21':'QUAD 3','P22':'QUAD 4','P23':'QUAD 3','P24':'QUAD 4',}
    
    comp['quad'] = comp['Well'].map(quad_384).fillna('missing')
       
    print(comp.head())
    st.table(comp.head())
    
    ####need a function here to strip out empty arrays from the data_stream - not a great idea
    #####stripping data out- function should perhaps sit in loop and check the mean of ROX - if less than 1000
    ####delete the whole array from the dataframe - might be a point during ArayaManger to run mean of ROX and add colution to 
    ####filter out and delete before changing index to odrder
        
    # these need to be wrapped in functions - currently outside functions for testing. Cache this and have it clear cache with
    #above clear all files button.
    
    #scoring - see if it is worth putting selection box for current cutoffs
    
    
    
  
    
    #comp.to_csv('file_uploader_check.csv')
    
   ####start producing plots - use st.plotly_chart() container to view in app - they can be downloaded to html as well - read docs on use######## 
    
    ####use this pands view to check changes to the datafrme are correct.
    
    
    #Display view of data for normalised N1N2 - order of processing over time with cutoff lines. 
   
    def fam_pro_qc(comp):
        figN1 = px.scatter(comp, x= 'order', y = 'norm_N_Cov' ,color = 'Result', title = 'N1 N2 Calls - normalised processinng view')

        figN1.add_trace(go.Scatter(
            y=[1, 1],
            x=[comp.order.min(), comp.order.max()],
            mode="lines+markers+text",
            name="Lower_1_Positive_Boundary",
            text=["1"],
            textposition="top center",
            line=dict(color="red", dash = 'dash')
            ))

        figN1.update_traces(marker_size=3)

        #figN1.add_trace(go.Scatter(
         #   y=[9, 9],
         #   x=[comp.order.min(), comp.order.max()],
         #   mode="lines+markers+text",
         #   name="Lower_9_Positive_Boundary",
         #   text=["9"],
          #  textposition="top center",
          #  line=dict(color="red", dash = 'dash')))

    

        figN1.update_xaxes(showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')
        figN1.update_yaxes(range=[0, 2.5],showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')
        st.plotly_chart(figN1, use_container_width=True)
    
   
    def RP_pro_QC(comp):
    
        
        fig1bbnbb = px.scatter(comp, x= 'order', y = 'norm_RNaseP', color = 'Result', title = 'normalised RNaseP processing view' )
        fig1bbnbb.update_traces(marker_size=3)
        fig1bbnbb.update_yaxes(range=[0, 1])
        fig1bbnbb.add_trace(go.Scatter(
            y=[0.2, 0.2],
            x=[comp.order.min(), comp.order.max()],
            mode="lines+markers+text",
            name="Lower_0.2_RP Detected_Boundary",
            text=["0.2"],
            textposition="top center",
            line=dict(color="red", dash = 'dash')))
        #fig1bbnbb.add_trace(go.Scatter(
        #    y=[2, 2],
        #    x=[comp.order.min(), comp.order.max()],
        #    mode="lines+markers+text",
        #    name="Lower_2_RP Detected_Boundary",
        #    text=["2"],
         #   textposition="top center",
         #   line=dict(color="blue", dash = 'dash')))
        st.plotly_chart(fig1bbnbb, use_container_width=True)
    
    
    
    
    
    def plot_roxCV(comp):
        ROX_mean = round(comp.ROX_RFU.mean())
        ROX_std = round(comp.ROX_RFU.std())
        CV = round(((ROX_std/ROX_mean)*100),1)
        CT = ("CV% "+ str(CV), "Mean "+ str(ROX_mean), "standard deviation "+str(ROX_std))
        print(CT)

        fig3a = px.scatter(comp, x= comp.order, y = comp.ROX_RFU ,color = comp.Result)
        fig3a.update_traces(mode='markers', marker_line_width=0.01, marker_size=2)
        fig3a.add_trace(go.Scatter(
            y=[ROX_mean, ROX_mean],
            x=[comp.order.min(), comp.order.max()],    
            mode="lines+text",
            name="Mean average",
            text=["Mean"],
            textposition="top center",
            line=dict(color="blue", dash = 'dash')))
        fig3a.add_trace(go.Scatter(
            y=[ROX_mean + (ROX_std * 3), ROX_mean + (ROX_std * 3)],
            x=[comp.order.min(), comp.order.max()],    mode="lines+markers+text",
            name="UCL - Upper +3 SD Cutoff",
            text=["UCL"],
            textposition="top center",
            line=dict(color="Red", dash = 'dash')))

        fig3a.add_trace(go.Scatter(
            y=[ROX_mean - (ROX_std * 3), ROX_mean - (ROX_std * 3)],
            x=[comp.order.min(), comp.order.max()],    
            mode="lines+text",
            name="-3SD ROX RFU LCL",
            text=["LCL"],
            textposition="top center",
            line=dict(color="Red", dash = 'dash')))
        fig3a.update_layout(title = 'ROX dispense processing ' + str(CT))
        fig3a.update_traces(marker_size=3)
        fig3a.update_yaxes(range=[4000, 12000])
        
        st.plotly_chart(fig3a, use_container_width=True)
        
    def plot_FAMCV(comp):
        ROX_mean = round(comp.FAM_RFU.mean())
        ROX_std = round(comp.FAM_RFU.std())
        CV = round(((ROX_std/ROX_mean)*100),1)
        CT = ("CV% "+ str(CV), "Mean "+ str(ROX_mean), "standard deviation "+str(ROX_std))
        print(CT)

        fig3a = px.scatter(comp, x= comp.order, y = comp.FAM_RFU ,color = comp.Result)
        fig3a.update_traces(mode='markers', marker_line_width=0.01, marker_size=2)
        fig3a.add_trace(go.Scatter(
            y=[ROX_mean, ROX_mean],
            x=[comp.order.min(), comp.order.max()],    
            mode="lines+text",
            name="Mean average",
            text=["Mean"],
            textposition="top center",
            line=dict(color="blue", dash = 'dash')))
        fig3a.add_trace(go.Scatter(
            y=[ROX_mean + (ROX_std * 3), ROX_mean + (ROX_std * 3)],
            x=[comp.order.min(), comp.order.max()],    mode="lines+markers+text",
            name="UCL - Upper +3 SD Cutoff",
            text=["UCL"],
            textposition="top center",
            line=dict(color="Red", dash = 'dash')))

        fig3a.add_trace(go.Scatter(
            y=[ROX_mean - (ROX_std * 3), ROX_mean - (ROX_std * 3)],
            x=[comp.order.min(), comp.order.max()],    
            mode="lines+text",
            name="-3SD ROX RFU LCL",
            text=["LCL"],
            textposition="top center",
            line=dict(color="Red", dash = 'dash')))
        fig3a.update_layout(title = 'ROX dispense processing ' + str(CT))
        fig3a.update_traces(marker_size=3)
        fig3a.update_yaxes(range=[4000, 12000])
        
        st.plotly_chart(fig3a, use_container_width=True) 
    
    #Display ROX vs FAM plot for over chemical performace vs dispense
    
    def roxfam(comp):
        ROX_mean = round(comp.ROX_RFU.mean())
        ROX_std = round(comp.ROX_RFU.std())
        CV = round(((ROX_std/ROX_mean)*100),1)
        CT = ("CV% "+ str(CV), "Mean "+ str(ROX_mean), "standard deviation "+str(ROX_std))
        print(CT)
        
        figroxfam = px.scatter(comp, x= 'ROX_RFU', y = 'FAM_RFU' ,color = 'Result', title = 'N1 N2 detection Performance vs dispense')
        figroxfam.update_traces(marker_size=3)
        figroxfam.update_xaxes(range=[4000, 10000])
        figroxfam.update_yaxes(range=[0, 20000])


        figroxfam.add_trace(go.Scatter(
            x=[ROX_mean - (ROX_std * 3), ROX_mean - (ROX_std * 3)],
            y=[20000, 0],
            mode="lines",
            name="RFU Lower Cutoff Limit",
            #text=["LCL"],
            text=["ROX -3SD lower cutoff"],
            textposition="top center",
            line=dict(color="Red", dash = 'dash')
            ))
        st.plotly_chart(figroxfam, use_container_width=True)
    
    def cluster(comp):
        fig2b = px.scatter(comp, x= 'norm_RNaseP', y = 'norm_N_Cov',color = 'Result', title = 'Cluster Processing view')
        fig2b.update_xaxes(range=[0, 1])
        fig2b.update_yaxes(range=[0, 2])
        fig2b.update_traces(marker_size=3)
        st.plotly_chart(fig2b, use_container_width=True)
    
    
    
    
    
    
    #ctrl_qc_table = testdf.groupby(['date_time','Result'])['norm_N_Cov','FAM_RFU', 'ROX].agg('mean', 'std')
    #ctrl_qc_table = testdf.groupby(['date_time','Result'])['norm_N_Cov','FAM_RFU', 'ROX].agg('mean', 'std')
    
    #print(ctrl_qc_table)   
    
    def ctrl_view(testdf):
        
        figdt = px.scatter(testdf, x='date_time', y='norm_N_Cov', color = 'Result')
        figdt.update_yaxes(range=[0, 2.5])
        figdt.update_traces(marker_size=3)
        figdt.add_trace(go.Scatter(
            y=[1.3, 1.3],
            x=[testdf.date_time.min(), testdf.date_time.max()],
            mode="lines+markers+text",
            name="Val mean",
            text=["Mean"],
            textposition="top center",
            line=dict(color="red", dash = 'dash')))
        figdt.add_trace(go.Scatter(
            y=[13.5, 13.5],
            x=[testdf.date_time.min(), testdf.date_time.max()],
            mode="lines+markers+text",
            name="3SD low",
            text=["-3SD"],
            textposition="top center",
            line=dict(color="yellow", dash = 'dash')))
        figdt.add_trace(go.Scatter(
            y=[15.2, 15.2],
            x=[testdf.date_time.min(), testdf.date_time.max()],
            mode="lines+markers+text",
            name="3SD high",
            text=["+3SD"],
            textposition="top center",
            line=dict(color="yellow", dash = 'dash')))
        st.plotly_chart(figdt, use_container_width=True)
    
    #percentiles
    def Q25(x):
        return x.quantile(0.25)

    def Q50(x):
        return x.quantile(0.5)

    def Q75(x):
        return x.quantile(0.75)
    
    
    #def valid_array(df):
     #   df = df.groupby('Run_ID')['ROX_RFU'].mean()
      #  st.table(df)
        
    #valid_array(comp)    
    
    
    def stats_nFAM(df):
        
        df['norm_N_Cov'] = df['norm_N_Cov'].abs()
        stats_nFAM = df.groupby(['Run_ID','quad', 'Result'])['norm_N_Cov'].agg(['count', 'mean','min', 'std', 'max']).fillna('-')
        
  
        
        CI95_hi_nFAM = []
        CI95_lo_nFAM = []
        CV_run_nFAM = []
  

        for i in stats_nFAM.index:
            c,m,s,t,v =(stats_nFAM.loc[i])
            CI95_hi_nFAM.append(m + 1.96*s/math.sqrt(c))
            CI95_lo_nFAM.append(m - 1.96*s/math.sqrt(c))
            CV_run_nFAM.append(100 - (s/m*100))
    
        stats_nFAM['CI 95% low nFAM'] = CI95_lo_nFAM
        stats_nFAM['CI 95_hi_nFAM'] = CI95_hi_nFAM
        stats_nFAM['CV%_nFAM'] = CV_run_nFAM
        #stats_nFAM['%Percent_detected'] = result['N1N2_detected'] / TOT*100
        return(stats_nFAM)
        
        
    def stats_nCFO(df):
        
        df['norm_RNAseP'] = df['norm_RNaseP'].abs()
        stats_nCFO = df.groupby(['Run_ID', 'quad', 'Result'])['norm_RNaseP'].agg(['count', 'mean','min', 'std', 'max']).fillna('-')
        
  
        
        CI95_hi_nCFO = []
        CI95_lo_nCFO = []
        CV_run_nCFO = []
  

        for i in stats_nCFO.index:
            c,m,s,t,v =(stats_nCFO.loc[i])
            CI95_hi_nCFO.append(m + 1.96*s/math.sqrt(c))
            CI95_lo_nCFO.append(m - 1.96*s/math.sqrt(c))
            CV_run_nCFO.append(100 - (s/m*100))
    
        stats_nCFO['CI 95% low nCFO'] = CI95_lo_nCFO
        stats_nCFO['CI 95_hi_nCFO'] = CI95_hi_nCFO
        stats_nCFO['CV%_nCFO'] = CV_run_nCFO
        #stats_nFAM['%Percent_detected'] = result['N1N2_detected'] / TOT*100
        return(stats_nCFO)
        
    
    
    #app layout - charts are already produced above - this allows arangment of charts in order to 
    #order the chart layout as required. Place any further charts /  table genertion above this line.
    #Place the st containers below this line to arrange them as required. leave tops headeders above all 
    #scripts so that the buttons etc.. all continue to work with the flow of Streamlit. 
    
    st.subheader('Nexar processing view - Trend Data')
    
    fam_pro_qc(comp)
    
    RP_pro_QC(comp)
    
    plot_roxCV(comp)
    
    st.subheader('Nexar processing - Clustering data view')
   

    col1, col2 = st.columns(2)

    with col1:
        roxfam(comp)
        
    with col2:
        cluster(comp)
        
    
    st.subheader('Accuplex Control View - QC control data')
    
   


   
    testdf = comp[(comp.control == 'A1500')]
    
    ctrl_view(testdf)
    
    nFAM_data = stats_nFAM(comp)
    nCFO_data = stats_nCFO(comp)
    
    st.dataframe(nFAM_data.astype(str))
    st.dataframe(nCFO_data.astype(str))
    
    
    def convert_df(df):
     # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv().encode('utf-8')

    all_data = convert_df(comp)

    st.download_button(
        label="Download all data as CSV",
         data=all_data,
         file_name='araya_all_data_ox.csv',
         mime='text/csv')
         

    
         
    nFAM = convert_df(nFAM_data)

    st.download_button(
        label="Download nFAM CSV",
         data=nFAM,
         file_name='araya_ox_out.csv',
         mime='text/csv')
    
    nCFO = convert_df(nCFO_data)

    st.download_button(
        label="Download nFAM CSV",
         data=nCFO,
         file_name='araya_ox_out.csv',
         mime='text/csv')
    
    
    
    #st.dataframe(round(ctrl_qc_table),0)
    
    
    
    

###function parser - parse araya files - instatiated with ArayaManager - functions can be accessed with . notation
class WellDataManager:
    """Parent class for importing data from any instrument"""
    def __init__(self,
                 files=None,
                 dfs=None,
                 run_name_split_loc=1,
                 date_time_split_loc=3,
                 group_name="",
                 replicate_name=""):

        # project attributes
        self.file_names = files
        self.df_list = dfs
        self.group_name = group_name
        self.replicate_name = replicate_name
        self.run_name = ""
        self.date_time = ""
        self.split_char_loc = run_name_split_loc
        self.run_df = pd.DataFrame()
        self.group_df = pd.DataFrame()

        # csv read attributes
        self.tables = 1
        self.index_column = [0]
        self.header_row = [1]
        self.row_count = [8]
        self.column_names = ['Row_ID', 'Col_ID', 'Value']

    def concatenate_dataframes(self, file_names):
        if file_names is not None:
            self.file_names = file_names
            self.build_dataframes(self.file_names)
            self.group_df = pd.concat([self.group_df, self.run_df], ignore_index=True)
        else:
            for each_file in self.file_names:
                self.file_names = file_names
                self.build_dataframes(each_file)
                self.group_df = pd.concat([self.group_df, self.run_df], ignore_index=True)
        return self.group_df

    def build_dataframes(self, each_file):
        self.read_file(each_file)
        self.coerce_numeric_values()
        self.run_df['Group_ID'] = self.group_name
        #self.run_df['File_root'] = each_file
        self.run_df['Run_ID'] = self.run_name
        self.run_df['date_time'] = self.date_time
        # print(self.run_df)

    def coerce_numeric_values(self):
        # may be used used to force values to numeric.  Not applicable to all instruments
        # may be used used to force values to numeric.  Not applicable to all instruments
        pass

    def read_file(self, file_name):
        """Reads Initial Data from CSV file"""
        df = pd.read_csv(file_name, header=self.header_row, nrows=self.row_count, index_col=self.index_column, dtype = 'Int32', thousands=',')
        df = df.stack()
        self.run_df = df.reset_index()
        self.run_df.columns = self.column_names

    def get_run_name(self, file_name):
        """Splits string to get run name from file name."""
        self.run_name = file_name[:self.split_char_loc]
    def get_date_time(self, file_name):
        self.date_time = file_name[-(self.split_char_loc+1):-3]
        #print(self.run_name)


class ArayaManager(WellDataManager):
    """Class that handles Well Data Data"""
    def __init__(self,
                 files=None,
                 dfs=None,
                 run_name_split_loc=6,
                 date_time_split_loc=1,
                 group_name="",
                 replicate_name="",
                 dyes=None,
                 separate_column=True):
        super().__init__(files,
                         dfs,
                         run_name_split_loc,
                         group_name,
                         replicate_name)

        if dyes is None:
            dyes = ['FAM', 'VIC', 'ROX']

        # Ayara-specific items
        self.separate_column_per_dye = separate_column
        self.channel_df = pd.DataFrame()
        self.dyes = dyes
        self.channels = ['CH1', 'CH2', 'CH3']
        self.length = 14

        # csv read attributes
        self.tables = 3
        self.index_column = ["<>", "<>", "<>"]
        self.header_row = [5, 23, 41]
        self.row_count = [16, 16, 16]

        if self.separate_column_per_dye:
            # TODO: generalize for other dye names
            self.column_names = ['Row_ID', 'Col_ID', 'FAM_RFU', 'VIC_RFU', 'ROX_RFU']
        else:
            self.column_names = ['Row_ID', 'Col_ID', 'RFU', 'Channel', 'Dye']

    def read_each_channel(self, df, ch):
        """Reads Individual Channel Data from CSV file"""

        # Need to shift to get rid of annoying '<>'.  Otherwise won't parse correctly.

        if '<>' in df:
            df = df.shift(periods=1, axis='columns')
            df.drop('<>', axis=1, inplace=True)

        # Stack df for various dyes and add additional columns
        df = df.stack()
        self.channel_df = df.reset_index()

        # For separate columns for each dye, rename RFU columns.  pd.concat() method does the rest!
        if self.separate_column_per_dye:
            self.channel_df.columns = self.column_names[0:3]
            self.channel_df.rename(columns={'FAM_RFU': f'{self.dyes[ch]}_RFU'},
                                   inplace=True)

        # case to stack all dyes into common RFU and Dye channels.
        else:
            self.channel_df['Channel'] = self.channels[ch]
            self.channel_df['Dye'] = self.dyes[ch]
            self.channel_df.columns = self.column_names


    def read_file(self, file_name):
        """Reads Each Channel Data from CSV file"""

        # loops through the 3 channel tables in the csv output files.
        self.run_df = pd.DataFrame()

        df = pd.read_csv(file_name,
                         header=self.header_row[0],
                         na_values="<>", thousands=',')
        dyes = ['FAM', 'VIC', 'ROX']
        dfs = {dye: pd.DataFrame() for dye in dyes}
        ranges = [(0, 16), (18, 34), (36, 52)]
        for i, dye in zip(ranges, dyes):
            dfs[dye] = df[i[0]:i[1]]

        for ch in range(self.tables):
            dye = dyes[ch]
            self.read_each_channel(dfs[dye], ch)

            # case to have separate columns for each dye
            if self.separate_column_per_dye:
                self.channel_df = self.channel_df[self.channel_df.columns.difference(self.run_df.columns)]
                self.run_df = pd.concat([self.run_df, self.channel_df], axis=1)

            # case to stack all dyes into common RFU and Dye channels.
            else:
                self.run_df = pd.concat([self.run_df, self.channel_df], ignore_index=True)

        # Force columns to correct order.  Fixed bug with concat of separate dye columns.
        self.run_df = self.run_df[self.column_names]

    def get_run_name(self, file_name):
        """Splits string to get run name from file name."""
        self.run_name = file_name[-(self.split_char_loc+4):-4]
        """Splits string of file name to collect date time"""
    def get_date_time(self, file_name):
        self.date_time = file_name[-(self.split_char_loc+50):-38]
        #add in statement to deal with appended num in duplicated upload files.
        if len(self.date_time) > self.length:
            self.date_time = file_name[-(self.split_char_loc+50):-40]
        
        


if __name__ == "__main__":
    main()
