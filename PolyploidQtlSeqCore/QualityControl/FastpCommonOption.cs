using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// Fastp共通オプション
    /// </summary>
    internal class FastpCommonOption
    {
        // 常に使用するFastpオプション
        private static readonly string _CutTailArg = "-3";
        private static readonly string _detectAdapterPEArg = "--detect_adapter_for_pe";
        private static readonly string _jsonReportArg = "-j /dev/null";


        /// <summary>
        /// Fastp共通オプションを作成する。
        /// </summary>
        /// <param name="readLengthRequired">リード最低長</param>
        /// <param name="nBaseLimit">N塩基リミット</param>
        /// <param name="baseQuality">塩基クオリティ</param>
        /// <param name="cutTailMeanQuality">3'末端平均クオリティ</param>
        /// <param name="cutTailWindowSize">3'末端トリムウインドウサイズ</param>
        /// <param name="outDir">出力ディレクトリ</param>
        /// <param name="threadNumber">スレッド数</param>
        public FastpCommonOption(ReadLengthRequired readLengthRequired, NBaseLimit nBaseLimit, BaseQuality baseQuality,
            CutTailMeanQuality cutTailMeanQuality, CutTailWindowSize cutTailWindowSize, OutputDirectory outDir, ThreadNumber threadNumber)
        {
            ReadLengthRequired = readLengthRequired;
            NBaseLimit = nBaseLimit;
            BaseQuality = baseQuality;
            CutTailMeanQuality = cutTailMeanQuality;
            CutTailWindowSize = cutTailWindowSize;
            OutputDirectory = outDir;
            ThreadNumber = threadNumber;
        }

        /// <summary>
        /// トリム後リード長の最低値
        /// </summary>
        public ReadLengthRequired ReadLengthRequired { get; }

        /// <summary>
        /// N塩基数の上限値
        /// </summary>
        public NBaseLimit NBaseLimit { get; }

        /// <summary>
        /// 塩基クオリティのしきい値
        /// </summary>
        public BaseQuality BaseQuality { get; }

        /// <summary>
        /// 3'末端トリム時の平均クオリティしきい値
        /// </summary>
        public CutTailMeanQuality CutTailMeanQuality { get; }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズ
        /// </summary>
        public CutTailWindowSize CutTailWindowSize { get; }

        /// <summary>
        /// 使用するスレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDirectory { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <param name="inputFilePair">入力Fastpペア</param>
        /// <returns>Fastp引数</returns>
        public string ToFastpArg(FastqFilePair inputFilePair)
        {
            var outputFastqPair = OutputDirectory.ToOutputFastqFilePair(inputFilePair);
            var htmlReportFile = new FastpHtmlReportFile(OutputDirectory, inputFilePair);

            var fastpArgs = new[]
            {
                $"-i {inputFilePair.Fastq1Path} -I {inputFilePair.Fastq2Path}",
                outputFastqPair.ToFastpArg(),
                ReadLengthRequired.ToFastpArg(),
                NBaseLimit.ToFastpArg(),
                BaseQuality.ToFastpArg(),
                _CutTailArg,
                CutTailWindowSize.ToFastpArg(),
                CutTailMeanQuality.ToFastpArg(),
                ThreadNumber.ToFastpArg(),
                _detectAdapterPEArg,
                htmlReportFile.ToFastpArg(),
                _jsonReportArg
            };

            return string.Join(" ", fastpArgs);
        }
    }
}
