using PolyploidQtlSeqCore.MappingAndVariantCall;

namespace PolyploidQtlSeq.Options.Pipeline
{
    internal class VariantCallPipelineSettingValueDTO : IVariantCallPipelineSettingValue
    {
        public VariantCallPipelineSettingValueDTO()
        {
            ReferenceSequence = "";
            ThreadNumber = 0;
            Parent1Dir = "";
            Parent2Dir = "";
            Bulk1Dir = "";
            Bulk2Dir = "";
            ChrSizeThreshold = 0;
            AnalysisChrNames = "";
            MinBq = 0;
            MinMq = 0;
            AdjustMq = 0;
            OutputDir = "";
            SnpEffMaxHeap = 0;
            SnpEffConfigFile = "";
            SnpEffDatabaseName = "";
        }


        public string ReferenceSequence { get; init; }

        public int ThreadNumber { get; init; }

        public string Parent1Dir { get; init; }

        public string Parent2Dir { get; init; }

        public string Bulk1Dir { get; init; }

        public string Bulk2Dir { get; init; }

        public int ChrSizeThreshold { get; init; }

        public string AnalysisChrNames { get; init; }

        public int MinMq { get; init; }

        public int MinBq { get; init; }

        public int AdjustMq { get; init; }

        public string OutputDir { get; init; }

        public int SnpEffMaxHeap { get; init; }

        public string SnpEffConfigFile { get; init; }

        public string SnpEffDatabaseName { get; init; }
    }
}
