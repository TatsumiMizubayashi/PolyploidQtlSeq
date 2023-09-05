using PolyploidQtlSeq.Options.QtlAnalysis;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// QTL-Seqパイプラインオプション値 インターフェース
    /// </summary>
    internal interface IQtlSeqPipelineOptionValue : IQtlSeqAnalysisOptionValue
    {
        public string ReferenceSequence { get; set; }

        public string Parent1Dir { get; set; }

        public string Parent2Dir { get; set; }

        public string Bulk1Dir { get; set; }

        public string Bulk2Dir { get; set; }

        public int ChrSizeThreshold { get; set; }

        public string AnalysisChrNames { get; set; }

        public int MinMq { get; set; }

        public int MinBq { get; set; }

        public int AdjustMq { get; set; }

        public int SnpEffMaxHeap { get; set; }

        public string SnpEffConfigFile { get; set; }

        public string SnpEffDatabaseName { get; set; }

        /// <summary>
        /// VariantCallPipelineSettingValueを作成する。
        /// </summary>
        /// <returns>VariantCallPipelineSettingValue</returns>
        internal VariantCallPipelineSettingValueDTO CreateVariantCallPipelineSettingValue()
        {
            return new VariantCallPipelineSettingValueDTO()
            {
                ReferenceSequence = ReferenceSequence,
                Parent1Dir = Parent1Dir,
                Parent2Dir = Parent2Dir,
                Bulk1Dir = Bulk1Dir,
                Bulk2Dir = Bulk2Dir,
                ChrSizeThreshold = ChrSizeThreshold,
                AnalysisChrNames = AnalysisChrNames,
                MinMq = MinMq,
                MinBq = MinBq,
                AdjustMq = AdjustMq,
                OutputDir = OutputDir,
                ThreadNumber = ThreadNumber,
                SnpEffMaxHeap = SnpEffMaxHeap,
                SnpEffConfigFile = SnpEffConfigFile,
                SnpEffDatabaseName = SnpEffDatabaseName
            };
        }
    }
}
