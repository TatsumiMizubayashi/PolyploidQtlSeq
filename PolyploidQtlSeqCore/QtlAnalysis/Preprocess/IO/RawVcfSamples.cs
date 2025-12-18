namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// RawVcfSampleコレクション
    /// </summary>
    internal class RawVcfSamples
    {
        private static readonly char[] _delimiter = [':'];
        private static readonly int _parent1Index = 0;
        private static readonly int _parent2Index = 1;
        private static readonly int _bulk1Index = 2;
        private static readonly int _bulk2Index = 3;

        private const string GT_KEY = "GT";
        private const string AD_KEY = "AD";


        /// <summary>
        /// RawVcfSampleコレクションインスタンスを作成する。
        /// </summary>
        /// <param name="altAlleles">Altアレル配列</param>
        /// <param name="format">FORMAT</param>
        /// <param name="sampleValues">サンプル値配列</param>
        public RawVcfSamples(string[] altAlleles, string format, string[] sampleValues)
        {
            if (altAlleles.Length == 0) throw new ArgumentException(null, nameof(altAlleles));
            if (string.IsNullOrWhiteSpace(format)) throw new ArgumentException(null, nameof(format));
            if (sampleValues.Length != 4) throw new ArgumentException(null, nameof(sampleValues));

            var (adIndex, gtIndex) = GetIndex(format);

            Parent1 = CreateRawSample(sampleValues[_parent1Index], adIndex, gtIndex);
            Parent2 = CreateRawSample(sampleValues[_parent2Index], adIndex, gtIndex);
            Bulk1 = CreateRawSample(sampleValues[_bulk1Index], adIndex, gtIndex);
            Bulk2 = CreateRawSample(sampleValues[_bulk2Index], adIndex, gtIndex);
            IsMultiAllele = altAlleles.Length > 1;
        }

        /// <summary>
        /// 親1を取得する。
        /// </summary>
        public RawVcfSample Parent1 { get; }

        /// <summary>
        /// 親2を取得する。
        /// </summary>
        public RawVcfSample Parent2 { get; }

        /// <summary>
        /// Bulk1を取得する。
        /// </summary>
        public RawVcfSample Bulk1 { get; }

        /// <summary>
        /// Bulk2を取得する。
        /// </summary>
        public RawVcfSample Bulk2 { get; }

        /// <summary>
        /// マルチアレルかどうかを取得する。
        /// </summary>
        public bool IsMultiAllele { get; }

        /// <summary>
        /// VcfParent1を取得する。
        /// </summary>
        /// <param name="allAlleles">全アレル</param>
        /// <returns>VcfParent1</returns>
        public VcfParent1 GetVcfParent1(string[] allAlleles)
        {
            var allele = Parent1.GT.ToAllele(allAlleles);

            return new VcfParent1(Parent1.GT, Parent1.AD, allele);
        }

        /// <summary>
        /// VcfParent2を取得する。
        /// </summary>
        /// <param name="allAlleles">全アレル</param>
        /// <returns>VcfParent2</returns>
        public VcfParent2 GetVcfParent2(string[] allAlleles)
        {
            var allele = Parent2.GT.ToAllele(allAlleles);

            return new VcfParent2(Parent2.GT, Parent2.AD, allele);
        }

        /// <summary>
        /// VcfBulk1を取得する。
        /// </summary>
        /// <param name="allAlleles">全アレル</param>
        /// <returns>VcfParent1</returns>
        public VcfBulk1 GetVcfBulk1(string[] allAlleles)
        {
            var allele = Bulk1.GT.ToAllele(allAlleles);

            return new VcfBulk1(Bulk1.GT, Bulk1.AD, allele);
        }

        /// <summary>
        /// VcfBulk2を取得する。
        /// </summary>
        /// <param name="allAlleles">全アレル</param>
        /// <returns>VcfParent1</returns>
        public VcfBulk2 GetVcfBulk2(string[] allAlleles)
        {
            var allele = Bulk2.GT.ToAllele(allAlleles);

            return new VcfBulk2(Bulk2.GT, Bulk2.AD, allele);
        }

        /// <summary>
        /// AD indexとGT indexを取得する。
        /// </summary>
        /// <param name="format">FORMAT</param>
        /// <returns>(ADindex, GTindex)</returns>
        /// <exception cref="ArgumentException">ADまたはGTが存在しない</exception>
        private static (int adIndex, int gtIndex) GetIndex(string format)
        {
            var formatItems = format.Split(_delimiter).ToList();
            var adIndex = formatItems.IndexOf(AD_KEY);
            if (adIndex == -1) throw new ArgumentException($"{AD_KEY} does not exist.");

            var gtIndex = formatItems.IndexOf(GT_KEY);
            if (gtIndex == -1) throw new ArgumentException($"{GT_KEY} does not exist.");

            return (adIndex, gtIndex);
        }

        /// <summary>
        /// RawVcfSampleインスタンスを作成する。
        /// </summary>
        /// <param name="value">サンプル値</param>
        /// <param name="adIndex">AD index</param>
        /// <param name="gtIndex">GT index</param>
        /// <returns>RawVcfSample</returns>
        private static RawVcfSample CreateRawSample(string value, int adIndex, int gtIndex)
        {
            var valueItems = value.Split(_delimiter);
            var ad = new AD(valueItems[adIndex]);
            var gt = new GT(valueItems[gtIndex]);

            return new RawVcfSample(ad, gt);
        }
    }
}
