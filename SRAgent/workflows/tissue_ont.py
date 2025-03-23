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
from SRAgent.agents.tissue_ont import create_tissue_ont_agent

# classes
class UBERON_IDS(BaseModel):
    #descriptions: List[str] = Field(description="Each unique tissue description provided in the input.")
    ids: List[str] = Field(description="The selected Uberon ID (UBERON:XXXXXXX)")

# functions
def create_tissue_ont_workflow(
    model_name: Optional[str]=None,
    return_tool: bool=True,
) -> Callable:
    # create model
    model = set_model(model_name=model_name, agent_name="tissue_ont")

    # set tools
    tools = [
        create_tissue_ont_agent(),
    ]
  
    # state modifier
    state_mod = "\n".join([
        "# Introduction",
        " - You are a helpful senior bioinformatician assisting a researcher with a task involving classifying one or more tissue.",
        " - You will be provided with a free text description of the tissues.",
        " - Your task is to categorize the tissues based on the Uberon ontology.",
        " - You must find the single most suitable Uberon ontology term that best describes the tissue description.",
        " - You have a set of tools that can help you with this task.",
        "# Tools",
        " - create_tissue_ont_agent: Use this tool to find the most suitable Uberon ontology term that best describes the tissue description.",
        "# Workflow",
        " 1. Identify each unique tissue description in the input.",
        "   - For example, 'brain cortex; eye lens; aortic valve;' should be split into the following separate descriptions:",
        "     - brain cortex",
        "     - eye lens",
        "     - aortic valve",
        " 2. For each description (e.g., \"brain cortex\"), use the create_tissue_ont_agent tool to find the most suitable Uberon ontology term.",
        "   - You MUST use the create_tissue_ont_agent tool for EACH tissue description."
        "# Response",
        " - Provide the list of most suitable Uberon ontology IDs (UBERON:XXXXXXX) that best describe each tissue description.",
    ])
    # create agent
    agent = create_react_agent(
        model=model,
        tools=tools,
        state_modifier=state_mod,
        response_format=UBERON_IDS,
    )

    # create tool
    @tool
    async def invoke_tissue_ont_workflow(
        messages: Annotated[List[BaseMessage], "Messages to send to the Tissue Ontology workflow"],
        config: RunnableConfig,
    ) -> Annotated[dict, "Response from the Tissue Ontology workflow"]:
        """
        Invoke the Tissue Ontology workflow with a message.
        The Tissue Ontology workflow will annotate each tissue description with the most suitable Uberon term.
        """
        # Invoke the agent with the message
        result = await agent.ainvoke({"messages" : messages}, config=config)
        return {"tissue_ontology_term_ids" : result["structured_response"].ids}
    return invoke_tissue_ont_workflow


# main 
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv(override=True)

    async def main():
        # create workflow
        workflow = create_tissue_ont_workflow()

        # Example 1: Complex tissue description example
        print("\n=== Example 1: Complex tissue description example ===")
        msg = "Categorize the following tissues: the thin layer of epithelial cells lining the alveoli in lungs; brain cortex; eye lens"
        input = {"messages": [HumanMessage(content=msg)]}
        ontology_term_ids = await workflow.ainvoke(input)
        print(f"Result for complex description: {','.join(ontology_term_ids['tissue_ontology_term_ids'])}")

    # run
    asyncio.run(main())