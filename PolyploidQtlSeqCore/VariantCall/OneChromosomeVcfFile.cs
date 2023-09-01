using PolyploidQtlSeqCore.QtlAnalysis.Chr;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// 染色体1本分のVCFファイル
    /// </summary>
    internal class OneChromosomeVcfFile
    {
        /// <summary>
        /// 染色体1本分のVCFファイルを作成する。
        /// </summary>
        /// <param name="filePath">VCFファイルのPath</param>
        /// <param name="chromosome">対象染色体</param>
        public OneChromosomeVcfFile(string filePath, Chromosome chromosome)
        {
            Chr = chromosome;
            Path = filePath;
        }

        /// <summary>
        /// 染色体情報を取得する。
        /// </summary>
        public Chromosome Chr { get; }

        /// <summary>
        /// VCFファイルPathを取得する。
        /// </summary>
        public string Path { get; }

        /// <summary>
        /// VCF indexファイルを作成する。
        /// </summary>
        /// <returns></returns>
        public async ValueTask CreateIndexFile()
        {
            await Tabix.RunAsync(Path);
        }

        /// <summary>
        /// VCFファイルとindexファイルを削除する。
        /// </summary>
        public void Delete()
        {
            if (File.Exists(Path)) File.Delete(Path);

            var indexFilePath = Path + ".tbi";
            if (File.Exists(indexFilePath)) File.Delete(indexFilePath);
        }
    }
}
