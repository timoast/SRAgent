# import
## batteries
import os
import asyncio
from typing import Annotated, List, Optional, Callable
## 3rd party
from pydantic import BaseModel, Field
from Bio import Entrez
from langchain_core.tools import tool
from langgraph.prebuilt import create_react_agent
from langchain_core.runnables.config import RunnableConfig
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
## package
from SRAgent.agents.utils import set_model
from SRAgent.tools.tissue_ont import query_vector_db, get_neighbors, query_uberon_ols

# classes
class UBERON_ID(BaseModel):
    id: str = Field(description="The Uberon ID (UBERON:XXXXXXX)")

# functions
def create_tissue_ont_agent(
    model_name: Optional[str]=None,
    return_tool: bool=True,
) -> Callable:
    # create model
    model = set_model(model_name=model_name, agent_name="tissue_ont")

    # set tools
    tools = [
        query_vector_db,
        get_neighbors,
        query_uberon_ols,
    ]
  
    # state modifier
    state_mod = "\n".join([
        "# Introduction",
        " - You are a helpful senior bioinformatician assisting a researcher with a task involving classifying tissues.",
        " - You will be provided with a free text description of a tissue.",
        " - Your task is to categorize the tissue based on the Uberon ontology.",
        " - You must find the single most suitable Uberon ontology term that best describes the tissue.",
        " - You have a set of tools that can help you with this task.",
        "# Tools",
        " - query_vector_db: Perform a semantic search on a vector database to find Uberon terms related to the target tissue. The database contains a collection of tissue descriptions and their corresponding Uberon terms.",
        " - get_neighbors: Get the neighbors of a given Uberon term in the Uberon ontology. Useful for finding adjacent terms in the ontology.",
        " - query_uberon_ols: Query the Ontology Lookup Service (OLS) for Uberon terms matching the search term.",
        "# Workflow",
        " Step 1: Start with the query_vector_db tools to find the most similar Uberon terms to the target tissue.",
        " Step 2: Use the get_neighbors tool on the Uberon terms returned in Step 1 to help find the most suitable term.",
        "   - ALWAYS use the get_neighbors tool to explore more the terms adjacent to the terms returned in Step 1.",
        " Step 3: Repeat steps 1 and 2 to find the most suitable Uberon term.",
        "   - ALWAYS perform between 1 and 3 iterations to find the most suitable term.",
        " Step 4: If you are uncertain about which term to select, use the query_uberon_ols tool to help find the most suitable term."
        "# Response",
        " - Provide the most suitable Uberon ontology ID (UBERON:XXXXXXX) that best describes the tissue.",
    ])
    # create agent
    agent = create_react_agent(
        model=model,
        tools=tools,
        state_modifier=state_mod,
        response_format=UBERON_ID,
    )

    # return agent instead of tool
    if not return_tool:
        return agent

    # create tool
    @tool
    async def invoke_tissue_ont_agent(
        messages: Annotated[List[BaseMessage], "Messages to send to the Tissue Ontology agent"],
        config: RunnableConfig,
    ) -> Annotated[dict, "Response from the Tissue Ontology agent"]:
        """
        Invoke the Tissue Ontology agent with a message.
        The Tissue Ontology agent will annotate a tissue description with the most suitable Uberon term.
        """
        # Invoke the agent with the message
        result = await agent.ainvoke({"messages" : messages}, config=config)
        return {"tissue_ontology_term_id" : result["structured_response"].id}
    return invoke_tissue_ont_agent

# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv(override=True)

    async def main():
        # create entrez agent
        agent = create_tissue_ont_agent(return_tool=False)
    
        # # Example 1: Simple tissue example
        # print("\n=== Example 1: Simple tissue example ===")
        # msg = "Categorize the following tissue: brain"
        # input = {"messages": [HumanMessage(content=msg)]}
        # result = await agent.ainvoke(input)
        # print(f"Result for 'brain': {result['structured_response'].id}")
        
        # Example 2: More specific tissue example
        # print("\n=== Example 2: More specific tissue example ===")
        # msg = "Categorize the following tissue: hippocampus"
        # input = {"messages": [HumanMessage(content=msg)]}
        # result = await agent.ainvoke(input)
        # print(f"Result for 'hippocampus': {result['structured_response'].id}")
        
        # # Example 3: Complex tissue description example
        print("\n=== Example 3: Complex tissue description example ===")
        msg = "Categorize the following tissue: the thin layer of epithelial cells lining the alveoli in lungs"
        input = {"messages": [HumanMessage(content=msg)]}
        result = await agent.ainvoke(input)
        print(f"Result for complex description: {result['structured_response'].id}")
        
    asyncio.run(main())