using PolyploidQtlSeqCore.QtlAnalysis.Chr;

namespace PolyploidQtlSeqCore.QtlAnalysis.Mapping
{
    /// <summary>
    /// 全サンプルのBAMファイル
    /// </summary>
    internal class AllSampleBamFiles
    {
        /// <summary>
        /// 全サンプルのBAMファイルを作成する。
        /// </summary>
        /// <param name="parent1Bam">Parent1のBAM</param>
        /// <param name="parent2Bam">Parent2のBAM</param>
        /// <param name="bulk1Bam">Bulk1のBAM</param>
        /// <param name="bulk2Bam">Bulk2のBAM</param>
        public AllSampleBamFiles(BamFile parent1Bam, BamFile parent2Bam, BamFile bulk1Bam, BamFile bulk2Bam)
        {
            Parent1BamFile = parent1Bam;
            Parent2BamFile = parent2Bam;
            Bulk1BamFile = bulk1Bam;
            Bulk2BamFile = bulk2Bam;
        }

        /// <summary>
        /// Parent1のBAMファイル
        /// </summary>
        public BamFile Parent1BamFile { get; }

        /// <summary>
        /// Parent2のBAMファイル
        /// </summary>
        public BamFile Parent2BamFile { get; }

        /// <summary>
        /// Bulk1のBAMファイル
        /// </summary>
        public BamFile Bulk1BamFile { get; }

        /// <summary>
        /// Bulk2のBAMファイル
        /// </summary>
        public BamFile Bulk2BamFile { get; }

        /// <summary>
        /// 解析対象染色体を取得する。
        /// </summary>
        /// <param name="option">解析対象染色体オプション</param>
        /// <returns>解析対象染色体</returns>
        public Chromosome[] GetAnalysisChrs(AnalysisChrOption option)
        {
            var bamHeader = Parent1BamFile.ToHeader();
            var allChrs = bamHeader.ToChromosomes();

            return allChrs.Get(option);
        }
    }
}
