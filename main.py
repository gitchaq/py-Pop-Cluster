import streamlit as st
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import MultiPolygon, MultiLineString,Polygon,Point,MultiPoint,LineString
import numpy as np
import folium
from rtree import index
from shapely.geometry import Point
from streamlit_folium import folium_static
from streamlit_option_menu import option_menu
from collections import Counter
from pathlib import Path
from shapely.geometry import shape
from streamlit_extras import faker, card, word_importances, app_logo
from PIL import Image
from pathlib import Path



# st.set_page_config(page_title="Areal Density",
#                    page_icon="https://aptaworks.com/wp-content/themes/aptaworks/images/favicon.ico",
#                    layout="wide")
app_logo.add_logo(
    "https://aptaworks.com/wp-content/themes/aptaworks/images/logo.png")
# dfg_pop=gpd.read_file('popdensity.geojson')
# dfg_market=gpd.read_file('market_clean.geojson')
hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 


selected2 = option_menu("Areal Density", ["Home", "Map", "About"],
                        icons=['house', 'cloud-upload', 'info-circle'], menu_icon="window-dock",
                        default_index=0, orientation="horizontal")
# selected
# st.title("Areal Density")



container = [st.container()]
container2 = [st.container()]
colArray1 = [5,2]
colArray2 = [5,5,5,5]
colArray3 = [5,5]
colArray4 = [5,10]

inp01, inp02, inp03,inp04 = container2[0].columns(colArray2)
input01, input02, input03,input04 = container2[0].columns(colArray2)
head1,head2=container2[0].columns(colArray4)
content0,content = container[0].columns(colArray1)

content01, content02, content03,content04 = container2[0].columns(colArray2)
overall,overall2 = container2[0].columns(colArray3)

#content0.title("")
dir_pop = Path(__file__).parent / 'population_fulldim_demo.geojson'
dir_market = Path(__file__).parent / 'market_fulldim_demo.geojson'
@st.cache(allow_output_mutation=True)
def load_data(dirpop,dirmarket):
    # with gzip.open(dir_pop, 'rb') as f:
    #     with tempfile.NamedTemporaryFile(mode='wb', delete=False) as tmp:
    #         tmp.write(f.read())
    #         tmp.flush()
    #         dfg_population = gpd.read_file(tmp.name, driver='GeoJSON')
    # #gpd.read_file("kontur_boundaries_ID_20220407.gpkg") 
    #with gzip.open(dirpop, 'rb') as f:
    dfg_population = gpd.read_file(dirpop).drop_duplicates(subset='geometry', keep='first').reset_index(drop=True)
    #dfg_population = dfg_population.drop_duplicates(subset='geometry', keep='first').reset_index(drop=True)
    #dfg_population=gpd.read_file(dir_pop).drop_duplicates(subset='geometry', keep='first').reset_index(drop=True)
    dfg_market=gpd.read_file(dirmarket).drop_duplicates(subset='geometry', keep='first').reset_index(drop=True)
    # Set up Rtree index
    return dfg_population,dfg_market
#dfg_clusters = gpd.GeoDataFrame(columns=['point_id', 'population', 'geometry'])

dfg_population_raw, dfg_market_raw = load_data(dir_pop,dir_market)
prov=dfg_population_raw['Province'].unique()



if selected2 == 'Home':
    # st_faker = faker.get_streamlit_faker(seed=42)
    st.subheader('Welcome')

    st.write(
        f'<hr style="background-color: {"#803df5"}; margin-top: 0;'
        ' margin-bottom: 0; height: 3px; border: none; border-radius: 3px;">',
        unsafe_allow_html=True,
    )
    text = """&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This application is designed to understand map of
                 |Population Density|. It uses a native hexagonal map data
                 from Indonesia Population Density and Indonesia Administration Area
                 and Traditional Market. We generate into more Intuitive and insightful
                 to understand better in Geospatial Analytics |This map only generate
                 Few administrative Location| additionally, We have the full Indonesia Data.
                 Overall, this app aims to generate |comprehensive geospatial analytics| for a specific case.
                 """
    words = text.split('|')
    desc = word_importances.format_word_importances(
        words,
        importances=(0, 0.5, 0, 0.5, 0, 0.5, 0),  # fmt: skip
    )
    st.write(desc, unsafe_allow_html=True)
    image = Image.open(Path(__file__).parent / 'density.png')
    st.image(image,use_column_width='always')


elif selected2 == 'Map':
    prov_select = head1.selectbox('Select a Province', prov,index=4)

    # if prov_select==prov:
    #     dfg_market=dfg_market_raw.loc[[dfg_market_raw['Province']==prov_select]].reset_index(drop=True)
    #     dfg_population=dfg_population_raw.loc[[dfg_population_raw['Province']==prov_select]].reset_index(drop=True)
    if prov_select is not None:
        dfg_market = dfg_market_raw.loc[dfg_market_raw['Province'].isin([prov_select])].reset_index(drop=True)
        dfg_population=dfg_population_raw.loc[dfg_population_raw['Province'].isin([prov_select])].reset_index(drop=True)
    else:
        dfg_market = dfg_market_raw.copy()
        dfg_population=dfg_population_raw.copy()

    idx = index.Index()
    for i, geom in enumerate(dfg_population.geometry):
        idx.insert(i, geom.bounds)

    #st.write(dfg_population.geometry.type.unique())
    #st.write(dfg_market.geometry.type.unique())
    MIN_POP = head2.slider('Minimum Population', 2000, 50000, 10000, 2000)

    max_pol = 200
    # Loop over all market points
    # used_pop=[]
    # used_geom=[]
    results = []
    # results_polygon=[]
    used_indexes = []
    for i, market_id in enumerate(dfg_market.index):
        #st.write(market_id)
        #st.write(i)
        # Find nearest polygon
        market_point = dfg_market.loc[market_id, 'geometry']
        #st.write(market_point)
        nearest_polygons = [j for j in idx.nearest(market_point.bounds, num_results=max_pol)]
        #st.write(nearest_polygons)
        nearest_polygons.sort(key=lambda j: market_point.distance(dfg_population.geometry.iloc[j]))
        #st.write(nearest_polygons)
        #nearest_polygon = nearest_polygons[0]
        # # Find nearest polygon
        # market_point = dfg_market.loc[market_id, 'geometry']
        # nearest_polygons = [j for j in idx.nearest(market_point.bounds, num_results=max_pol)]
        # nearest_polygons.sort(key=lambda j: market_point.distance(dfg_population.geometry.iloc[j]))
        # if not nearest_polygons:
        #     continue
        # nearest_polygon = nearest_polygons[0]
        #st.write(nearest_polygons)  
        # Calculate total population within threshold distance
        total_pop = 0
        total_area=0
        pop_indexes = []
        for j, pop_poly in enumerate(dfg_population.geometry.iloc[nearest_polygons]):
            total_pop += dfg_population.loc[nearest_polygons[j], 'Population']
            total_area+=450
            pop_indexes.append(nearest_polygons[j])
            if total_pop >= MIN_POP:
                break
        # Add polygons and total population to results
        result = {
            'id': market_id+1,
            'geometry': dfg_population.iloc[pop_indexes].unary_union,
            'population': total_pop,
            'area(m2)':total_area
        }
        #used_index.append(pop_indexes)
        #st.write(used_index)
        geo_polygon=[shapely.wkt.loads(str(geom)) for geom in dfg_population.iloc[pop_indexes].geometry]
        pop_polygon=[pop for pop in dfg_population.iloc[pop_indexes].Population]
        result_polygon = {
            'geometry': geo_polygon,
            'population': pop_polygon
            #'area(m2)':450
        }
        results.append(result)
        used_indexes.extend(pop_indexes)
    dfg_result = gpd.GeoDataFrame(results)
    dfg_result = dfg_result.set_crs(dfg_population.crs)
    dfg_result=dfg_result.to_crs("EPSG:4326")


    # unique_indexes =set(used_indexes)
    duplicate_indexes = [index for index, count in Counter(used_indexes).items() if count > 1]
    uniqused_indexes=list(set(used_indexes))
    used_mask = dfg_population.index.isin(uniqused_indexes)
    #used_mask = dfg_population.index.isin(used_indexes)
    duplicated_pop=dfg_population.iloc[duplicate_indexes]
    #uncovered_pop= dfg_population.iloc[~used_indexes]
    uncovered_pop = dfg_population.loc[~used_mask]
    uncovered_pop=uncovered_pop.drop_duplicates(subset='geometry', keep='first')
    # new_dictg=[]
    # new_dictp=[]
    # for item in results_polygon:
    #     newgeom = item["geometry"]
    #     newpop = item["population"]
    #     new_dictg.concat(newgeom)
    #     new_dictp.concat(newpop)
    # new_dict={
    #         'geometry': list(new_dictg),
    #         'population': new_dictp
    #         #'area(m2)':450
    #     }
    # results_polygon_all = new_dict
    # st.write(results_polygon)    
    #st.write(results_polygon) 
    # result = {
    #     'geometry': results_polygon_all,
    #     'population': 450,
    #     'area(m2)':450
    # }        
    # results_polygon_all = [[elem[0] if isinstance(elem, tuple) else elem for elem in sublist] for sublist in results_polygon]
    #results_polygon=[(geom) for geom, index in results_polygon]
    # Create result GeoDataFrame
    # st.write(results_polygon_all)
    #st.write(tes)
    #st.write(results)

    # for geometry_list in results_polygon_all['geometry']:
    #     for i, geometry in enumerate(geometry_list):
    #         results_polygon_all['geometry'][geometry_list][i] = geometry[0]

    # for population_list in results_polygon_all['population']:
    #     for i, population in enumerate(population_list):
    #         results_polygon_all['population'][population_list][i] = population[0]
    # results_polygon_all=gpd.GeoDataFrame(results_polygon_all,geometry='geometry').reset_index(drop=True)
    # st.write(results_polygon_all)
    # duplicated_geometry = results_polygon_all[results_polygon_all.duplicated(subset='geometry', keep=False)].copy()
    # duplicated_geometry=gpd.GeoDataFrame(duplicated_geometry)


    # Print result
    # st.write(dfg_result)
    # Set up Rtree index


    # st.write(duplicated_pop)
    # st.write(dfg_market)

    uncovered_pop=uncovered_pop.loc[~uncovered_pop['geometry'].isna()].reset_index(drop=True)

    idx = index.Index()
    for i, geom in enumerate(uncovered_pop.geometry):
        idx.insert(i, geom.bounds)
    used_ucindexes=[]
    # Loop over all uncovered polygons

    uncovered_results = []
    for i, uncovered_id in enumerate(uncovered_pop.index):
        # Check if polygon has already been used
        if uncovered_id in used_ucindexes:
            continue
        
        # Find polygons within threshold distance
        uncovered_poly = uncovered_pop.loc[uncovered_id, 'geometry']
        nearest_polygons = [j for j in idx.nearest(uncovered_poly.bounds, num_results=max_pol)]
        nearest_polygons.sort(key=lambda j: uncovered_poly.distance(uncovered_pop.geometry.iloc[j]))

        # Calculate total population within threshold distance
        total_pop = 0
        pop_indexes = []
        total_area=0
        for j, pop_poly in enumerate(uncovered_pop.geometry.iloc[nearest_polygons]):
            # Check if polygon has already been used
            if nearest_polygons[j] in used_ucindexes:
                continue
            total_area +=450
            total_pop += uncovered_pop.loc[nearest_polygons[j], 'Population']
            pop_indexes.append(nearest_polygons[j])
            if total_pop >= MIN_POP:
                break

        # Add polygons and total population to results
        result = {
            'id': len(uncovered_results)+1,
            'geometry': uncovered_pop.iloc[pop_indexes].unary_union,
            'population': total_pop,
            'area(m2)':total_area
        }
        uncovered_results.append(result)
        used_ucindexes.extend(pop_indexes)

    # Create result GeoDataFrame
    dfg_uncovered = gpd.GeoDataFrame(uncovered_results)
    dfg_uncovered = dfg_uncovered.set_crs(dfg_population.crs)
    dfg_uncovered = dfg_uncovered.to_crs("EPSG:4326")
    dfg_uncovered=dfg_uncovered.loc[~dfg_uncovered['geometry'].isna()]

    style1 = {'fillColor': '#82231b', 'color': '#82231b', 'lineColor': '#82231b'}
    style2 = {'fillColor': '#e2652d', 'color': '#e2652d', 'lineColor': '#e2652d'}
    map_ = folium.Map(location=[dfg_result.centroid.y.mean(), dfg_result.centroid.x.mean()], zoom_start=10)
    folium.GeoJson(dfg_uncovered.to_json(),
                style_function=lambda x:style2,
                #fill_opacity=1,
                tooltip=folium.GeoJsonTooltip(fields=['id', 'population','area(m2)'])).add_to(map_)
    folium.GeoJson(dfg_result.to_json(),
                tooltip=folium.GeoJsonTooltip(fields=['id', 'population','area(m2)'])).add_to(map_)
    folium.GeoJson(duplicated_pop.to_json(),
                style_function=lambda x:style1,
                #fill_opacity=1,
                tooltip=folium.GeoJsonTooltip(fields=['id', 'Population'])).add_to(map_)
    folium.GeoJson(dfg_market.to_json(),
                tooltip=folium.GeoJsonTooltip(fields=['Province', 'name','Kelompok Komoditas Utama','Klasifikasi','Perkiraan Jumlah Pedagang','Waktu Operasi'])).add_to(map_)

    folium_static(map_,width=1000,height=700)

    # #uncovered_pop
    # folium.GeoJson(uncovered_pop.to_json(),
    #                style_function=lambda x:style2,
    #                #fill_opacity=1,
    #                tooltip=folium.GeoJsonTooltip(fields=['id', 'Population'])).add_to(map_)
    # folium_static(map_,width=1500)
    # #uncovered_pop
    # folium.GeoJson(dfg_population.to_json(),
    #                style_function=lambda x:style2,
    #                #fill_opacity=1,
    #                tooltip=folium.GeoJsonTooltip(fields=['id', 'Population'])).add_to(map_)
    # folium_static(map_,width=1500)
    st.write('Red Polygon --> Intersect area between Market')
    st.write('Blue Polygon --> Supplied Market Area')
    st.write('Orange Polygon --> Potential Uncovered Area')

elif selected2 == 'About':
    st_faker = faker.get_streamlit_faker(seed=42)
    st_faker.subheader()
    st_faker.markdown()
